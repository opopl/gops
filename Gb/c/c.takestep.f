C   GMIN: A program for finding global minima!{{{
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of GMIN.
C
C   GMIN is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   GMIN is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C!}}}

      SUBROUTINE TAKESTEP(NP)

      USE COMMONS

      IMPLICIT NONE

      DOUBLE PRECISION DPRAND, RANDOM, XMASS, YMASS, ZMASS, LOCALSTEP, DUMMY2, CDIST(NATOMS), RDOTN, XL, YL, ZL,
     1                 DIST(3*NATOMS), DMAX, VMAX, VMIN, VMAX2, EXPDECAY(NATOMS), XC, YC, ZC, ANGLE, COST, SINT,
     2                 THETA, PHI, PI, DUMMY, CMDIST(NATOMS), CMMAX, RANDOMX, RANDOMY, RANDOMZ, RVEC(3), TX, TY, TZ,
     &                 DELX, DELY, DELZ, SLENGTH, RPROJ, X(3*NATOMS), DUMMYGRAD(3*NATOMS), DUMMYE, THETA2, FRPISQ, OPOTEL
      PARAMETER (PI=3.141592654D0)
      INTEGER J1, J2, JMAX, NP, J3, JMAX2, RANATOM, INDEXD(NATOMS), NCOOPDONE, ISTAT, NDUMMY, NTRIES, NTRIESMAX

      
C  Calling CENTRE if NORESET is .TRUE. can lead to problems with COORDSO containing an atom
C  outside the permitted radius. Then it may be impossible to take a step that keeps all the
C  atoms inside.
C
      FRPISQ = 4.D0*PI*PI
      NTRIESMAX=100

C     IF (CENT.AND.(.NOT.SEEDT)) CALL CENTRE2(COORDS(1:3*NATOMS,NP)) ! COORDS might have been shifted by symmetry
      IF ((.NOT.NORESET).AND.(.NOT.PERMOPT).AND.(.NOT.DIFFRACTT).AND.(.NOT.BLNT).AND.(.NOT.PERIODIC)
     &     .AND.(.NOT.GAUSST).AND.(.NOT.(CSMT.AND.(.NOT.SYMMETRIZECSM))).AND.(.NOT.PERCOLATET)) THEN
C
C	csw34> CHECK NOTHING HAS MOVED OUTSIDE THE CONTAINER RADIUS 
C
         DO J1=1,NATOMS
            IF ((.NOT.RIGID).OR.(J1.LE.NATOMS/2)) THEN
               J2=3*J1
               DUMMY2=COORDS(J2-2,NP)**2+COORDS(J2-1,NP)**2+COORDS(J2,NP)**2
               IF (DUMMY2.GT.RADIUS) THEN
                  IF (AMBERT) THEN ! jmc49 We don't really want a container at all in amber9, but this bit of code is being used 
                                   ! to warn about drift in, for instance, the amber-md.
                     WRITE(MYUNIT,'(A,I5,5F20.10)')'J1,RAD,D2,x,y,z=',J1,RADIUS,DUMMY2,COORDS(J2-2,NP),COORDS(J2-1,NP),COORDS(J2,NP)
                     WRITE(MYUNIT,'(A,I6)') 'WARNING initial coordinate outside container for particle ',J1
                     WRITE(MYUNIT,'(A)') 'NOT repositioning particle for AMBER9'
                     CYCLE ! the 1, natoms do loop.
                  END IF
                  WRITE(MYUNIT,'(A,I5,5F20.10)') 'J1,RAD,D2,x,y,z=',J1,RADIUS,DUMMY2,COORDS(J2-2,NP),COORDS(J2-1,NP),COORDS(J2,NP)
                  WRITE(MYUNIT,'(A,I6)') 'WARNING initial coordinate outside container - reposition particle ',J1
                  IF (DEBUG) THEN
                     WRITE(77,*) NATOMS
                     WRITE(77,*) ' '
                     WRITE(77,'(A3,3F20.10)') ('LA ',COORDS(3*(J3-1)+1,NP),COORDS(3*(J3-1)+2,NP),COORDS(3*(J3-1)+3,NP),J3=1,NATOMS)
                  ENDIF
C
C  Put it back in at the opposite end of a diameter
C
                  COORDS(J2-2,NP)=-COORDS(J2-2,NP)*0.8D0
                  COORDS(J2-1,NP)=-COORDS(J2-1,NP)*0.8D0
                  COORDS(J2,NP)=-COORDS(J2,NP)*0.8D0

               ENDIF
            ENDIF
         ENDDO
C
C	csw34> END OF RADIUS CHECK
C
      ENDIF
!
! COORDSO and VATO should now contain the correct COORDS and VAT even if
! these entries are changed by SYMMETRIZE or NEWRESET.
!
!     IF (.NOT.AMBERT) THEN
!        IF (ABS(COORDSO(1,NP)-COORDS(1,NP)).GT.1.0D-3) THEN
!           WRITE(MYUNIT,'(A,2G20.10)'),'takestep> WARNING - coordso will be changed: ',COORDSO(1,NP),COORDS(1,NP)
!        ENDIF
!        DO J1=1,3*(NATOMS-NSEED)
!           COORDSO(J1,NP)=COORDS(J1,NP)
!        ENDDO
!        DO J1=1,NATOMS
!           VATO(J1,NP)=VAT(J1,NP)
!        ENDDO
!     ENDIF

      IF (WENZEL) THEN
11       RANDOM=(DPRAND()-0.5D0)*2.0D0
         COORDS(1,NP)=COORDSO(1,NP)+STEP(NP)*RANDOM
         IF ((COORDS(1,NP).GT.1.0D0).OR.(COORDS(1,NP).LT.0.0D0)) GOTO 11
12       RANDOM=(DPRAND()-0.5D0)*2.0D0
         COORDS(2,NP)=COORDSO(2,NP)+STEP(NP)*RANDOM
         IF ((COORDS(2,NP).GT.1.0D0).OR.(COORDS(2,NP).LT.0.0D0)) GOTO 12
         RETURN
      ELSE IF (PERMOPT.OR.(CSMT.AND.(.NOT.SYMMETRIZECSM))) THEN
         RANDOM=(DPRAND()-0.5D0)*2.0D0
         ANGLE=RANDOM*STEP(NP) ! in radians
         COST=COS(ANGLE)
         SINT=SIN(ANGLE)
!
! Random rotation about the x axis.
!
         DO J1=1,NATOMS
            TY= COST*COORDS(3*(J1-1)+2,NP)+SINT*COORDS(3*(J1-1)+3,NP)
            TZ=-SINT*COORDS(3*(J1-1)+2,NP)+COST*COORDS(3*(J1-1)+3,NP)
            COORDS(3*(J1-1)+2,NP)=TY
            COORDS(3*(J1-1)+3,NP)=TZ
         ENDDO
         RANDOM=(DPRAND()-0.5D0)*2.0D0
         ANGLE=RANDOM*STEP(NP) ! in radians
         COST=COS(ANGLE)
         SINT=SIN(ANGLE)
!
! Random rotation about the y axis.
!
         DO J1=1,NATOMS
            TX= COST*COORDS(3*(J1-1)+1,NP)+SINT*COORDS(3*(J1-1)+3,NP)
            TZ=-SINT*COORDS(3*(J1-1)+1,NP)+COST*COORDS(3*(J1-1)+3,NP)
            COORDS(3*(J1-1)+1,NP)=TX
            COORDS(3*(J1-1)+3,NP)=TZ
         ENDDO
!
! Random rotation about the z axis.
!
         RANDOM=(DPRAND()-0.5D0)*2.0D0
         ANGLE=RANDOM*STEP(NP) ! in radians
         COST=COS(ANGLE)
         SINT=SIN(ANGLE)
         DO J1=1,NATOMS
            TX= COST*COORDS(3*(J1-1)+1,NP)+SINT*COORDS(3*(J1-1)+2,NP)
            TY=-SINT*COORDS(3*(J1-1)+1,NP)+COST*COORDS(3*(J1-1)+2,NP)
            COORDS(3*(J1-1)+1,NP)=TX
            COORDS(3*(J1-1)+2,NP)=TY
         ENDDO
!        WRITE(MYUNIT,'(A,I6)') 'takestep> coordinates rotated randomly as a rigid body for replica ',NP
!        IF (DEBUG) THEN
!           OPEN(UNIT=765,FILE='after_takestep.xyz',STATUS='UNKNOWN')
!           WRITE(765,'(I6)') NATOMS
!           WRITE(765,'(A,I6)') ' '
!           WRITE(765,'(A,3G20.10)') ('LA ',COORDS(3*(J2-1)+1,NP),COORDS(3*(J2-1)+2,NP),COORDS(3*(J2-1)+3,NP),J2=1,NATOMS)
!           CLOSE(765)
!        ENDIF

         RETURN
      ENDIF
C
C	csw34> FIND CENTRE OF COORDINATES (should be made a function!)
C
      XMASS=0.0D0; YMASS=0.0D0; ZMASS=0.0D0
      DO J1=1,NATOMS
         XMASS=XMASS+COORDS(3*(J1-1)+1,NP)
         YMASS=YMASS+COORDS(3*(J1-1)+2,NP)
         ZMASS=ZMASS+COORDS(3*(J1-1)+3,NP)
      ENDDO
      XMASS=XMASS/NATOMS; YMASS=YMASS/NATOMS; ZMASS=ZMASS/NATOMS
C
C  Find the most weakly bound atom, JMAX, the second most weakly bound atom, JMAX2,
C  and the pair energy of the most tightly bound atom, VMIN. An angular step is
C  taken for JMAX if its pair energy is > ASTEP*VMIN putting the atom at a radius of
C  DMAX (or CMMAX from CM of the cluster).
C

      DMAX=-1.0D0
      VMAX=-1.0D6
      VMAX2=-1.0D6
      VMIN=1.0D6
      CMMAX=-1.0D0
      DO J1=1,NATOMS
         J2=3*J1
         DIST(J1)= DSQRT( COORDS(J2-2,NP)**2+        COORDS(J2-1,NP)**2+        COORDS(J2,NP)**2)
         CMDIST(J1)=SQRT((COORDS(J2-2,NP)-XMASS)**2+(COORDS(J2-1,NP)-YMASS)**2+(COORDS(J2,NP)-ZMASS)**2)
!        IF (CMDIST(J1).GT.CMMAX) CMMAX=CMDIST(J1)
         IF ((CMDIST(J1).GT.CMMAX).AND.(J1.LE.NATOMS-NCORE(NP))) CMMAX=CMDIST(J1)
         IF (DIST(J1).GT.DMAX) DMAX=DIST(J1)
!        IF ((.NOT.SHELLMOVES(NP)).OR.(SHELLMOVES(NP).AND.(J1.LE.NATOMS-NCORE(NP)))) THEN
            IF (VAT(J1,NP).GT.VMAX) THEN
               VMAX=VAT(J1,NP)
               JMAX=J1
            ELSE IF ((VAT(J1,NP).LT.VMAX).AND.(VAT(J1,NP).GT.VMAX2)) THEN
               VMAX2=VAT(J1,NP)
               JMAX2=J1
            ENDIF
!        ENDIF
         IF (VAT(J1,NP).LT.VMIN) VMIN=VAT(J1,NP)
      ENDDO
C
C  If DECAY is true then select an atom at random, move this one randomly
C  by the maximum allowed amount, and move the others in random directions
C  of decaying magnitude, depending on how far they are from the chosen atom.
C
  
      IF (DECAY) THEN
9        RANATOM=NINT(0.5D0+NATOMS*DPRAND())
         IF (RANATOM.EQ.JMAX) GOTO 9 ! don't choose the atom that might undergo a surface move
         WRITE(MYUNIT,'(A,I6)') 'atom undergoing maximum displacement is number ',RANATOM
         DO J1=1,NATOMS-NSEED
            DUMMY=((COORDS(3*J1-2,NP)-COORDS(3*RANATOM-2,NP))**2+
     1             (COORDS(3*J1-1,NP)-COORDS(3*RANATOM-1,NP))**2+
     2             (COORDS(3*J1,NP)-  COORDS(3*RANATOM,NP))**2)
            EXPDECAY(J1)=EXP(-DECAYPARAM*DUMMY)
         ENDDO
      ENDIF
C
C  If MOVESHELL is true then we try a random angular move for all the frozen atoms 
C  about the centre of coordinates of the frozen set with probability SHELLPROB.
C
!     IF (SHELLMOVES(NP)) THEN
      IF (NCORE(NP).GT.0) THEN
         NSURFMOVES(NP)=NSURFMOVES(NP)+1
         XC=0.0D0; YC=0.0D0; ZC=0.0D0
         DO J1=NATOMS,NATOMS-NCORE(NP)+1,-1
            XC=XC+COORDS(3*(J1-1)+1,NP)
            YC=YC+COORDS(3*(J1-1)+2,NP)
            ZC=ZC+COORDS(3*(J1-1)+3,NP)
         ENDDO
         XC=XC/NCORE(NP); YC=YC/NCORE(NP); ZC=ZC/NCORE(NP)
         IF (DEBUG) WRITE(MYUNIT,'(A,3F12.4)') 'takestep> centre of coordinates for frozen atoms: ',XC, YC, ZC
         IF (MOVESHELLT.AND.(DPRAND().GT.(1.0D0-SHELLPROB))) THEN 
            WRITE(MYUNIT,'(A,I1,A,I8)') '[',NP,']takestep> shell move number ',NSURFMOVES(NP)

            RVEC(1)=(DPRAND()-0.5D0)*2.0D0
            RVEC(2)=(DPRAND()-0.5D0)*2.0D0
            RVEC(3)=(DPRAND()-0.5D0)*2.0D0
            DUMMY=SQRT(RVEC(1)**2+RVEC(2)**2+RVEC(3)**2)
            RVEC(1)=RVEC(1)/DUMMY; RVEC(2)=RVEC(2)/DUMMY; RVEC(3)=RVEC(3)/DUMMY
            ANGLE=DPRAND()*PI*2.0D0
            COST=COS(ANGLE)
            SINT=SIN(ANGLE)

            IF (DEBUG) WRITE(MYUNIT,'(A,I1,A,F10.2,A,3F12.4)') '[',NP,']takestep> angle=',ANGLE,' axis: ',RVEC(1:3)
C
C  Rotate all the non-core atoms through ANGLE about RVEC. Use rotation formula
C  from Goldstein p. 165.
C  
C           DO J1=NCORE(NP)+1,NATOMS
            DO J1=1,NATOMS-NCORE(NP)
               XL=COORDS(3*(J1-1)+1,NP); YL=COORDS(3*(J1-1)+2,NP); ZL=COORDS(3*(J1-1)+3,NP)
               DUMMY=SQRT((XL-XC)**2+(YL-YC)**2+(ZL-ZC)**2)
   
               RDOTN=(XL-XC)*RVEC(1)+(YL-YC)*RVEC(2)+(ZL-ZC)*RVEC(3)
               TX=(XL-XC)*COST + RVEC(1)*RDOTN*(1.0D0-COST)-((YL-YC)*RVEC(3)-(ZL-ZC)*RVEC(2))*SINT
               TY=(YL-YC)*COST + RVEC(2)*RDOTN*(1.0D0-COST)-((ZL-ZC)*RVEC(1)-(XL-XC)*RVEC(3))*SINT
               TZ=(ZL-ZC)*COST + RVEC(3)*RDOTN*(1.0D0-COST)-((XL-XC)*RVEC(2)-(YL-YC)*RVEC(1))*SINT
               IF (DUMMY.GT.0.1D0) THEN
                  COORDS(3*(J1-1)+1,NP)=(XC+TX)*(DUMMY+1.0D0)/DUMMY
                  COORDS(3*(J1-1)+2,NP)=(YC+TY)*(DUMMY+1.0D0)/DUMMY
                  COORDS(3*(J1-1)+3,NP)=(ZC+TZ)*(DUMMY+1.0D0)/DUMMY
               ELSE 
                  COORDS(3*(J1-1)+1,NP)=(XC+TX)
                  COORDS(3*(J1-1)+2,NP)=(YC+TY)
                  COORDS(3*(J1-1)+3,NP)=(ZC+TZ)
               ENDIF
            ENDDO
            IF (NSURFMOVES(NP).GE.SHELLMOVEMAX) THEN
               SHELLMOVES(NP)=.FALSE.
               NCORE(NP)=0
            ENDIF
            RETURN
         ENDIF
C        IF (NSURFMOVES(NP).GE.SHELLMOVEMAX) THEN
C           NSURFMOVES(NP)=-1
C           SHELLMOVES(NP)=.FALSE.
C           NCORE(NP)=0
C        ENDIF
333      CONTINUE
      ENDIF

C     IF (NSURFMOVES(NP).LT.0) THEN
C        NSURFMOVES(NP)=NSURFMOVES(NP)-1
C        IF (NSURFMOVES(NP).LE.-SHELLMOVEMAX) THEN
C           NSURFMOVES(NP)=0 
C        ENDIF
C     ENDIF
C
C  If COOP is true then select an atom at random, move this one randomly
C  by the maximum allowed amount, and move its NCOOP nearest neighbours by
C  the same displacement.
C  Note that the core atoms are ordered LAST, not FIRST !!
C
!     IF (COOP.AND.(NATOMS-NCORE(NP).GE.MAX(2,NCOOP))) THEN
      IF (COOP.AND.(NCORE(NP).GE.MAX(2,NCOOP))) THEN
!     IF (COOP) THEN
8        IF (NCORE(NP).GT.0) THEN
!           RANATOM=NINT(0.5D0+(NATOMS-NCORE(NP))*DPRAND())
!           RANATOM=NINT(0.5D0+(NATOMS-NSEED)*DPRAND())
            RANATOM=NATOMS-NCORE(NP)+NINT(0.5D0+NCORE(NP)*DPRAND())
         ELSE
            RANATOM=NINT(0.5D0+(NATOMS-NSEED)*DPRAND())
         ENDIF
         IF (RANATOM.EQ.JMAX) GOTO 8 ! don't choose the atom that might undergo a surface move
         IF (DEBUG) WRITE(MYUNIT,'(A,I6)') 'takestep> randomly selected atom for coop move is number ',RANATOM
         DO J1=1,NATOMS-NSEED

            CDIST(J1)=((COORDS(3*J1-2,NP)-COORDS(3*RANATOM-2,NP))**2+
     1                 (COORDS(3*J1-1,NP)-COORDS(3*RANATOM-1,NP))**2+
     2                 (COORDS(3*J1,NP)-  COORDS(3*RANATOM,NP))**2)
C           IF (SHELLMOVES(NP).AND.(J1.LE.NCORE(NP))) CDIST(J1)=1.0D100
C           IF (SHELLMOVES(NP).AND.(J1.GT.NATOMS-NCORE(NP))) CDIST(J1)=1.0D100
            IF ((NCORE(NP).GT.0).AND.(J1.GT.NATOMS-NCORE(NP))) CDIST(J1)=1.0D100
            INDEXD(J1)=J1
         ENDDO
         CALL SORT4(NATOMS-NSEED,NATOMS,CDIST,INDEXD)
!        PRINT *,'INDEXD=',INDEXD(1:NATOMS)
!        PRINT *,'CDIST=',CDIST(1:NATOMS)
!        CALL FLUSH(6,ISTAT)
         RANDOMX=(DPRAND()-0.5D0)*2.0D0
         RANDOMY=(DPRAND()-0.5D0)*2.0D0
         RANDOMZ=(DPRAND()-0.5D0)*2.0D0
         DUMMY2=SQRT(RANDOMX**2+RANDOMY**2+RANDOMZ**2)
         RANDOMX=RANDOMX/DUMMY2
         RANDOMY=RANDOMY/DUMMY2
         RANDOMZ=RANDOMZ/DUMMY2
         NCOOPDONE=0
      ENDIF
C
C	csw34> For the below loop, J1 is the atom number, and J2 is set to point to 3*J1, meaning
C	the z-coordinate of that atom. The y coordinate is then J2-1 and x is J2-2. 
C
      DO J1=1,NATOMS-NSEED
         IF (FROZEN(J1)) CYCLE
C        IF (SHELLMOVES(NP).AND.(J1.GT.NATOMS-NCORE(NP))) CYCLE
C        IF ((NCORE(NP).GT.0).AND.(J1.GT.NATOMS-NCORE(NP))) CYCLE
C        IF (NMOVE.GT.0) THEN
C           IF (1.0D0*NMOVE/1.0D0*NATOMS.LT.DPRAND()) GOTO 13
C        ENDIF
         IF ((COREFRAC.EQ.0.0D0).AND.(J1.GT.NATOMS-NCORE(NP))) CYCLE ! no point taking a zero step
         NTRIES=0
10       J2=3*J1
         IF ((NCORE(NP).GT.0).AND.(J1.GT.NATOMS-NCORE(NP))) THEN
            LOCALSTEP=STEP(NP)*COREFRAC ! smaller step for core
         ELSE
            LOCALSTEP=STEP(NP)
         ENDIF
         IF (RIGID.AND.(J1.GT.NATOMS/2)) THEN
            LOCALSTEP=0.0D0
            IF (OMOVE(NP)) LOCALSTEP=OSTEP(NP)
         ELSE IF (RIGID.AND.(J1.LE.NATOMS/2)) THEN
            LOCALSTEP=0.0D0
            IF (TMOVE(NP)) LOCALSTEP=STEP(NP)
         ENDIF
C
C  Angular move block.
C  If NORESET is .TRUE. then VAT won;t be set, so we should skip this block.
C
!        IF (J1.EQ.JMAX) WRITE(MYUNIT,'(A,I6,4F15.5)') 'JMAX,VAT,ASTEP(NP),VMIN,prod=',JMAX,VAT(J1,NP), 
!    &                                    ASTEP(NP),VMIN,ASTEP(NP)*VMIN
         IF (((VAT(J1,NP).GT.ASTEP(NP)*VMIN).AND.(J1.EQ.JMAX)).AND.(.NOT.RIGID).AND.(.NOT.BLNT).AND. 
     &         (.NOT.DIFFRACTT).AND.(.NOT.GAUSST).AND.(.NOT.PERCOLATET) 
     &        .AND.(.NOT.NORESET).AND.(.NOT.PERIODIC).AND.(.NOT.THOMSONT)
     &        .AND.(.NOT.((NCORE(NP).GT.0).AND.(J1.GT.NATOMS-NCORE(NP))))) THEN

            IF (DEBUG) WRITE(MYUNIT,'(A,I4,A,F12.4,A,F12.4,A,I4,A,F12.4)') 'angular move for atom ',J1, 
     &           ' V=',VMAX,' Vmin=',VMIN,' next most weakly bound atom is ',JMAX2,' V=',VMAX2

           THETA=DPRAND()*PI
           PHI=DPRAND()*PI*2.0D0
C
C  Evaporation is judged from the origin, not the centre of mass. We don't want the
C  angular move to cause evaporation. Obviously this will cause problems if we have a cluster that drifts
C  away from the origin up to the container radius.  
C
!          IF (SHELLMOVES(NP)) THEN ! different origin - significantly worse
!             COORDS(J2-2,NP)=XC+(CMMAX+1.0D0)*DSIN(THETA)*DCOS(PHI)
!             COORDS(J2-1,NP)=YC+(CMMAX+1.0D0)*DSIN(THETA)*DSIN(PHI)
!             COORDS(J2,NP)=  ZC+(CMMAX+1.0D0)*DCOS(THETA)
!          ELSE
              COORDS(J2-2,NP)=XMASS+(CMMAX+1.0D0)*DSIN(THETA)*DCOS(PHI)
              COORDS(J2-1,NP)=YMASS+(CMMAX+1.0D0)*DSIN(THETA)*DSIN(PHI)
              COORDS(J2,NP)=  ZMASS+(CMMAX+1.0D0)*DCOS(THETA)
!          ENDIF
           DUMMY=COORDS(J2-2,NP)**2+COORDS(J2-1,NP)**2+COORDS(J2,NP)**2
           IF (DUMMY.GT.RADIUS) THEN
              DUMMY=SQRT(RADIUS*0.99D0/DUMMY)
              COORDS(J2-2,NP)=COORDS(J2-2,NP)*DUMMY
              COORDS(J2-1,NP)=COORDS(J2-1,NP)*DUMMY
              COORDS(J2,NP)=COORDS(J2,NP)*DUMMY
           ENDIF
C
C  Angular move for the next most weakly bound atom doesn;t seem to help
C
C          ELSE IF ((.NOT.ROUND).AND.(((VAT(J1,NP).GT.ASTEP(NP)*VMIN).AND.(J1.EQ.JMAX2)).AND.(DZTEST))) THEN
C            IF (DEBUG) WRITE(MYUNIT,'(A,I4,A,F12.4,A,F12.4,A,I4,A,F12.4)') 
C     1                'angular move for atom ',J1,' V=',VMAX,' Vmin=',VMIN,' this is the next most weakly bound atom'
C            THETA=DPRAND()*PI
C            PHI=DPRAND()*PI*2.0D0
C
C            COORDS(J2-2,NP)=XMASS+(CMMAX+1.0D0)*DSIN(THETA)*DCOS(PHI)
C            COORDS(J2-1,NP)=YMASS+(CMMAX+1.0D0)*DSIN(THETA)*DSIN(PHI)
C            COORDS(J2,NP)=  ZMASS+(CMMAX+1.0D0)*DCOS(THETA)
C
C  End of angular move block.
C
         ELSE IF ((NATOMS-NSEED.EQ.1).AND.(NATOMS.GT.1)) THEN 
           IF (DEBUG) WRITE(MYUNIT,'(A,I4,A,F12.4,A,2F12.4)') 
     1                'angular move for atom ',J1,' V=',VAT(J1,NP),' Vmin, Vmax=',VMIN,VMAX
           THETA=DPRAND()*PI
           PHI=DPRAND()*PI*2.0D0
           COORDS(J2-2,NP)=XMASS+(CMMAX+1.0D0)*DSIN(THETA)*DCOS(PHI)
           COORDS(J2-1,NP)=YMASS+(CMMAX+1.0D0)*DSIN(THETA)*DSIN(PHI)
           COORDS(J2,NP)=  ZMASS+(CMMAX+1.0D0)*DCOS(THETA)
           DUMMY=COORDS(J2-2,NP)**2+COORDS(J2-1,NP)**2+COORDS(J2,NP)**2
           IF (DUMMY.GT.RADIUS) THEN
              DUMMY=SQRT(RADIUS*0.99D0/DUMMY)
              COORDS(J2-2,NP)=COORDS(J2-2,NP)*DUMMY
              COORDS(J2-1,NP)=COORDS(J2-1,NP)*DUMMY
              COORDS(J2,NP)=COORDS(J2,NP)*DUMMY
           ENDIF
         ELSE IF (DECAY) THEN
            RANDOMX=(DPRAND()-0.5D0)*2.0D0
            RANDOMY=(DPRAND()-0.5D0)*2.0D0
            RANDOMZ=(DPRAND()-0.5D0)*2.0D0
            DUMMY2=SQRT(RANDOMX**2+RANDOMY**2+RANDOMZ**2)
            COORDS(J2-2,NP)=COORDS(J2-2,NP)+LOCALSTEP*RANDOMX*EXPDECAY(J1)/DUMMY2
            COORDS(J2-1,NP)=COORDS(J2-1,NP)+LOCALSTEP*RANDOMY*EXPDECAY(J1)/DUMMY2
            COORDS(J2,NP)=  COORDS(J2,NP)+  LOCALSTEP*RANDOMZ*EXPDECAY(J1)/DUMMY2
!        ELSE IF (COOP.AND.(NATOMS-NCORE(NP).GE.MAX(2,NCOOP))) THEN
         ELSE IF (COOP.AND.(NCORE(NP).GE.MAX(2,NCOOP))) THEN
!        ELSE IF (COOP) THEN
            J2=3*INDEXD(J1)
            DUMMY=CDIST(J1)
            IF (NCORE(NP).GT.0) THEN
               J2=3*INDEXD(J1)
               DUMMY=CDIST(J1)
            ENDIF
!           WRITE(MYUNIT,'(A,2I6,F15.5,I6)') 'takestep> J1,J2/3,dist,NCOOPDONE=',J1,J2/3,DUMMY,NCOOPDONE
            IF ((NCOOPDONE.LE.NCOOP+1).AND.(DUMMY.LT.COOPCUT)) THEN
               NCOOPDONE=NCOOPDONE+1
               WRITE(MYUNIT,'(A,I6,A,F12.4)') 'takestep> coop move for atom ',J2/3,' cdist=',DUMMY
               COORDS(J2-2,NP)=COORDS(J2-2,NP)+LOCALSTEP*RANDOMX
               COORDS(J2-1,NP)=COORDS(J2-1,NP)+LOCALSTEP*RANDOMY
               COORDS(J2,NP)=  COORDS(J2,NP)+  LOCALSTEP*RANDOMZ
!
! Then a random move as well!
!
!              COORDS(J2-2,NP)=COORDS(J2-2,NP)+LOCALSTEP*(DPRAND()-0.5D0)*2.0D0
!              COORDS(J2-1,NP)=COORDS(J2-1,NP)+LOCALSTEP*(DPRAND()-0.5D0)*2.0D0
!              COORDS(J2,NP)=  COORDS(J2,NP)+  LOCALSTEP*(DPRAND()-0.5D0)*2.0D0
            ELSE
               COORDS(J2-2,NP)=COORDS(J2-2,NP)+LOCALSTEP*(DPRAND()-0.5D0)*2.0D0
               COORDS(J2-1,NP)=COORDS(J2-1,NP)+LOCALSTEP*(DPRAND()-0.5D0)*2.0D0
               COORDS(J2,NP)=  COORDS(J2,NP)+  LOCALSTEP*(DPRAND()-0.5D0)*2.0D0
            ENDIF
         ELSE IF (.NOT.FIXD) THEN
C
C  Tried changing to scale steps according to distance from CM. Maximum
C  allowed shift is linear with this distance. Worse.
C  Now try moving atoms according to how strongly bound they are.
C
           RANDOM=DPRAND()
C          WRITE(MYUNIT,'(A,3F20.10)') 'EFAC,RANDOM,EXP=',EFAC,RANDOM,EXP(-EFAC*(VAT(J1,NP)-VMAX)/(VMIN-VMAX))
C          PRINT*,'VMIN,VMAX,EFAC=',VMIN,VMAX,EFAC
           IF ((VMIN-VMAX.EQ.0.0D0).OR.(EFAC.EQ.0.0D0)) THEN
C             IF (.FALSE..AND.SHELLMOVES(NP)) THEN ! project out radial component and rescale to same length
                 IF (NCORE(NP).GT.0) THEN ! project out radial component and rescale to same length
                                          ! Only for shell atoms, not core, if they are moved at all.
                    DELX=(DPRAND()-0.5D0)*2.0D0*LOCALSTEP
                    DELY=(DPRAND()-0.5D0)*2.0D0*LOCALSTEP
                    DELZ=(DPRAND()-0.5D0)*2.0D0*LOCALSTEP
                    IF (J1.LE.NATOMS-NCORE(NP)) THEN
                       SLENGTH=SQRT(DELX**2+DELY**2+DELZ**2)
                       DUMMY=SQRT((COORDS(J2-2,NP)-XC)**2+(COORDS(J2-1,NP)-YC)**2+(COORDS(J2,NP)-ZC)**2)
                       RPROJ=(COORDS(J2-2,NP)-XC)*DELX+(COORDS(J2-1,NP)-YC)*DELY+(COORDS(J2,NP)-ZC)*DELZ
                       IF (DUMMY.NE.0.0D0) THEN
                          RPROJ=RPROJ/DUMMY
                          DELX=DELX-RPROJ*(COORDS(J2-2,NP)-XC)/DUMMY
                          DELY=DELY-RPROJ*(COORDS(J2-1,NP)-YC)/DUMMY
                          DELZ=DELZ-RPROJ*(COORDS(J2,NP)-ZC)/DUMMY
                          DUMMY=SQRT(DELX**2+DELY**2+DELZ**2)
                          IF (DUMMY.NE.0.0D0) THEN
                             DELX=DELX*SLENGTH/DUMMY
                             DELY=DELY*SLENGTH/DUMMY
                             DELZ=DELZ*SLENGTH/DUMMY
                          ENDIF
                       ENDIF
                    ENDIF
                    COORDS(J2-2,NP)=COORDS(J2-2,NP)+DELX
                    COORDS(J2-1,NP)=COORDS(J2-1,NP)+DELY
                    COORDS(J2,NP)=COORDS(J2,NP)+DELZ
                 ELSE
C
C	csw34> MAKE BASIC RANDOM CARTESIAN DISPLACEMENT MOVE
C              I've added the below IF line so that when MOVABLEATOMS is specified, cartesian displacements
C              specified by STEP will only be applied to atoms in that the 'movableatoms' file. It could be
C              that this needs changing into a new keyword to allow you to rotate bits and apply cartesian
C              moves to the whole system - something to keep in mind
C
C                if(.not.allocated(movableatomlist)) ALLOCATE(movableatomlist(1))
C                if(.not.allocated(movableatomlistlogical)) ALLOCATE(movableatomlistlogical(1))
C                IF (MOVABLEATOMST.AND.(.NOT.MOVABLEATOMLISTLOGICAL(J1))) CYCLE 
                 RANDOM=(DPRAND()-0.5D0)*2.0D0
                 IF (GAUSST) RANDOM=RANDOM/GKSMALL(J2-2) ! scale gauss steps by 1/GKSMALL
                 COORDS(J2-2,NP)=COORDS(J2-2,NP)+LOCALSTEP*RANDOM
                 RANDOM=(DPRAND()-0.5D0)*2.0D0
                 IF (GAUSST) RANDOM=RANDOM/GKSMALL(J2-1)
                 COORDS(J2-1,NP)=COORDS(J2-1,NP)+LOCALSTEP*RANDOM
                 RANDOM=(DPRAND()-0.5D0)*2.0D0
                 IF (GAUSST) RANDOM=RANDOM/GKSMALL(J2)
                 COORDS(J2,NP)=COORDS(J2,NP)+LOCALSTEP*RANDOM
C
C	csw34> END OF BASIC RANDOM CARTESIAN DISPLACEMENT MOVE
C
              ENDIF
           ELSE 
              RANDOM=(DPRAND()-0.5D0)*2.0D0
              DUMMY=1.0D0+EAMP*TANH(-2.0D0*EFAC*(VAT(J1,NP)-(VMAX+VMIN)/2.0D0)/(VMAX-VMIN))
C             COORDS(J2-2,NP)=COORDS(J2-2,NP)+LOCALSTEP*RANDOM*CMDIST(J1)/CMMAX
              COORDS(J2-2,NP)=COORDS(J2-2,NP)+LOCALSTEP*RANDOM*DUMMY
              RANDOM=(DPRAND()-0.5D0)*2.0D0
C             COORDS(J2-1,NP)=COORDS(J2-1,NP)+LOCALSTEP*RANDOM*CMDIST(J1)/CMMAX
              COORDS(J2-1,NP)=COORDS(J2-1,NP)+LOCALSTEP*RANDOM*DUMMY
              RANDOM=(DPRAND()-0.5D0)*2.0D0
C             COORDS(J2,NP)=COORDS(J2,NP)+LOCALSTEP*RANDOM*CMDIST(J1)/CMMAX
              COORDS(J2,NP)=COORDS(J2,NP)+LOCALSTEP*RANDOM*DUMMY
           ENDIF
C
C Stop atoms leaving the container in this step
C
           IF ((.NOT.PERIODIC).AND.(.NOT.AMBERT).AND.(.NOT.(RIGID.AND.((J1.GT.NATOMS/2)))).AND.(.NOT.BLNT)
     1     .AND.(.NOT.PERCOLATET).AND.(.NOT.DIFFRACTT).AND.(.NOT.THOMSONT).AND.(.NOT.GAUSST)) THEN
C          IF ((.NOT.PERIODIC).AND.(.NOT.AMBER).AND.(.NOT.(RIGID.AND.(LOCALSTEP.EQ.0.0D0))).AND.(.NOT.BLNT)) THEN
              DUMMY=COORDS(J2-2,NP)**2+COORDS(J2-1,NP)**2+COORDS(J2,NP)**2
C
C  Simply rescaling the radius of an atom that leaves the container will bias the sampling
C  of configuration space. However, we are not using takestep for bspt thermodynamic sampling!
C  So, put the atom back in the container on the other side!
C
              IF (DUMMY.GT.RADIUS) THEN
                 COORDS(J2-2,NP)=(SQRT(RADIUS)-0.5D0)*COORDS(J2-2,NP)/SQRT(DUMMY)
                 COORDS(J2-1,NP)=(SQRT(RADIUS)-0.5D0)*COORDS(J2-1,NP)/SQRT(DUMMY)
                 COORDS(J2,NP)=(SQRT(RADIUS)-0.5D0)*COORDS(J2,NP)/SQRT(DUMMY)
                 NTRIES=NTRIES+1
                 IF (NTRIES.GT.NTRIESMAX) THEN
                    WRITE(MYUNIT,'(A,I6,A)') 'takestep> WARNING atom ',J1,' persistently outside container'
                 ELSE
                    IF (DEBUG) WRITE(MYUNIT,'(A,I6,A,I6)') 'takestep> WARNING atom ',J1,' outside container, NTRIES=',NTRIES
                    GOTO 10
                 ENDIF
              ENDIF
           ENDIF
         ENDIF
         IF (TOSI.OR.WELCH) THEN
            DO J3=1,J1-1
               DUMMY=(COORDS(J2-2,NP)-COORDS(3*(J3-1)+1,NP))**2
     1              +(COORDS(J2-1,NP)-COORDS(3*(J3-1)+2,NP))**2
     2              +(COORDS(J2  ,NP)-COORDS(3*(J3-1)+3,NP))**2
               IF (DUMMY.LT.1.0D0) GOTO 10 
            ENDDO
         ENDIF
C13       CONTINUE
      ENDDO

      IF (FIXD) CALL HSMOVE(COORDS(1:3*NATOMS,1:NPAR),NP,NHSMOVE)
C
C  Preserve centre of mass if required.
C
      IF (CENT.AND.(.NOT.SEEDT)) CALL CENTRE2(COORDS(1:3*NATOMS,NP))
      IF (FIXCOM.AND.(.NOT.SEEDT)) CALL CENTRECOM(COORDS(:3*NATOMS,NP))

C     PRINT*,'NSYMREM in takestep=',NSYMREM
      IF (NSYMREM.GT.0) THEN
         CALL KEEPSYM(NP)
C        OPEN(UNIT=77,FILE='coords.latest.xyz',STATUS='OLD',POSITION='APPEND')
C        WRITE(77,'(I6)') NATOMS
C        WRITE(77,'(A)') ' '
C        WRITE(77,'(A2,2X,3G20.10)') ('LA',COORDS(3*(J2-1)+1,NP),COORDS(3*(J2-1)+2,NP),COORDS(3*(J2-1)+3,NP),J2=1,NATOMS)
C        CLOSE(77)
      ENDIF
!     DO J1=1,NATOMS
!        DO J2=J1+1,NATOMS
!           DUMMY=SQRT((COORDS(3*(J1-1)+1,NP)-COORDS(3*(J2-1)+1,NP))**2+
!    &              (COORDS(3*(J1-1)+2,NP)-COORDS(3*(J2-1)+2,NP))**2+
!    &              (COORDS(3*(J1-1)+3,NP)-COORDS(3*(J2-1)+3,NP))**2)
!           IF (DUMMY.LT.0.1D0) WRITE(MYUNIT,'(A,I6,A,I6,A,G20.10)') 'takestepC> WARNING, distance between atoms ',J1,' and ',
!    &                       J2,' is only ',DUMMY
!        ENDDO
!     ENDDO

!     DC430 >  
!     |gd351> added patchy 
      IF (DBPT.OR.DBPTDT.OR.DMBLMT.OR.LWOTPT.OR.MSSTOCKT.OR.MSTBINT.OR.NCAPT.OR.NPAHT.OR.NTIPT.OR.MULTPAHAT.OR.PAHAT.OR.STOCKAAT
     &    .OR.PAHW99T.OR.TDHDT.OR.PATCHY.OR.PAPT.OR.SILANET) THEN

         DO J1 = NATOMS/2 + 1, NATOMS
            J2      = 3*J1
            THETA2  = DOT_PRODUCT(COORDS(J2-2:J2,NP),COORDS(J2-2:J2,NP))
            IF (THETA2 > FRPISQ) THEN
               THETA2   = DSQRT(THETA2)
               THETA    = THETA2 - INT(THETA2/(2.D0*PI))*2.D0*PI
               COORDS(J2-2:J2,NP) = COORDS(J2-2:J2,NP)/THETA2*THETA
            ENDIF
         ENDDO

      ENDIF 

      RETURN
      END
C
      SUBROUTINE KEEPSYM(NP)
      USE commons
      IMPLICIT NONE
      INTEGER NP, J2, SYMINDEX(NATOMS), J3, J4
      DOUBLE PRECISION LCOORDS(3*NATOMS), NEWQ(3*NATOMS), SYMDELTA(3*NATOMS), DELTA(3*NATOMS),  SYMOP1(3,3)
      DOUBLE PRECISION STEPLENGTH, NEWSTEPLENGTH, ODIST1, ODIST2, DMIN, DUMMY
      LOGICAL ASSIGNED(NATOMS), BAD

      LCOORDS(1:3*NATOMS)=COORDSO(1:3*NATOMS,NP)
      DELTA(1:3*NATOMS)=COORDS(1:3*NATOMS,NP)-COORDSO(1:3*NATOMS,NP)
C     STEPLENGTH=SUM(DELTA(1:3*NATOMS)**2)
      SYMDELTA(1:3*NATOMS)=DELTA(1:3*NATOMS)

!     OPEN(UNIT=789,FILE='symstep.xyz',STATUS='UNKNOWN')
!     WRITE(789,'(I6)') NATOMS
!     WRITE(789,'(A)') 'coordso before symmetrising step'
!     WRITE(789,'(A3,3G20.10)') ('LA ',COORDSO(3*(J2-1)+1,NP),COORDSO(3*(J2-1)+2,NP),
!    &                                 COORDSO(3*(J2-1)+3,NP),J2=1,NATOMS-NCORE(NP))
!     WRITE(789,'(A3,3G20.10)') ('LB ',COORDSO(3*(J2-1)+1,NP),COORDSO(3*(J2-1)+2,NP),
!    &                                 COORDSO(3*(J2-1)+3,NP),J2=NATOMS-NCORE(NP)+1,NATOMS)
!     WRITE(789,'(I6)') NATOMS
!     WRITE(789,'(A)') 'coords before symmetrising step'
!     WRITE(789,'(A3,3G20.10)') ('LA ',COORDS(3*(J2-1)+1,NP),COORDS(3*(J2-1)+2,NP),
!    &                                 COORDS(3*(J2-1)+3,NP),J2=1,NATOMS-NCORE(NP))
!     WRITE(789,'(A3,3G20.10)') ('LA ',COORDS(3*(J2-1)+1,NP),COORDS(3*(J2-1)+2,NP),
!    &                                 COORDS(3*(J2-1)+3,NP),J2=NATOMS-NCORE(NP)+1,NATOMS)
C
C  New algorithm - choose the closest unclaimed atom in each case, so that
C  no tolerances are involved. 
C
      DO J2=1,NSYMREM
         BAD=.FALSE.
         SYMOP1(1:3,1:3)=SYMREM(J2,1:3,1:3)
         CALL MATMULV(NEWQ,LCOORDS,SYMOP1,NATOMS,3,3)
         ASSIGNED(1:NATOMS)=.FALSE.
         SYMINDEX(1:NATOMS)=0
         DO J3=1,NATOMS
            DMIN=1.0D100
            DO J4=1,NATOMS
               DUMMY=(LCOORDS(3*(J4-1)+1)-NEWQ(3*(J3-1)+1))**2+
     &               (LCOORDS(3*(J4-1)+2)-NEWQ(3*(J3-1)+2))**2+
     &               (LCOORDS(3*(J4-1)+3)-NEWQ(3*(J3-1)+3))**2
C              PRINT '(A,2I5,2G15.5)','J3,J4,DUMMY,DMIN=',J3,J4,DUMMY,DMIN
               IF (DUMMY.LT.DMIN) THEN
                  IF (ASSIGNED(J4)) THEN
C                    WRITE(MYUNIT,'(2(A,I5),A,F12.3)') 'WARNING closest atom ',J4,' to image atom ',J3,
C    &                                       ' already assigned, dist=',SQRT(DUMMY)
                  ELSE
                     IF (SYMINDEX(J3).GT.0) ASSIGNED(SYMINDEX(J3))=.FALSE.
                     SYMINDEX(J3)=J4
                     ASSIGNED(J4)=.TRUE.
                     DMIN=DUMMY
                  ENDIF
               ENDIF
            ENDDO
            IF (DEBUG.AND.(SQRT(DMIN).GT.SYMTOL5)) THEN
               WRITE(MYUNIT, '(2(A,I5),A,F12.3)') 'WARNING closest image to atom ',J3,' is ',SYMINDEX(J3), 
     &                                         ' distance=',SQRT(DMIN)
               BAD=.TRUE.
            ENDIF
         ENDDO
!        IF (BAD) THEN
!           OPEN(UNIT=1,FILE='keepsym.xyz',STATUS='UNKNOWN')
!           WRITE(1,*) NATOMS
!           WRITE(1,'(A)') 'LCOORDS'
!           DO J4=1,NATOMS
!              WRITE(1,'(A2,3X,3F20.10)') 'LA',LCOORDS(3*(J4-1)+1),LCOORDS(3*(J4-1)+2),LCOORDS(3*(J4-1)+3)
!           ENDDO
!           WRITE(1,*) NATOMS
!           WRITE(1,'(A)') 'NEWQ'
!           DO J4=1,NATOMS
!              WRITE(1,'(A2,3X,3F20.10)') 'LA',NEWQ(3*(J4-1)+1),NEWQ(3*(J4-1)+2),NEWQ(3*(J4-1)+3)
!           ENDDO
!           CLOSE(1)
!        ENDIF

         CALL MATMULV(NEWQ,DELTA,SYMOP1,NATOMS,3,3)
C
C  For each symmetry operation we average over the steps for the atoms
C  that are mapped on to after acting with the symmetry operation on the
C  original step vector.
C
         DO J3=1,NATOMS
            SYMDELTA(3*(J3-1)+1)=SYMDELTA(3*(J3-1)+1)+NEWQ(3*(SYMINDEX(J3)-1)+1)
            SYMDELTA(3*(J3-1)+2)=SYMDELTA(3*(J3-1)+2)+NEWQ(3*(SYMINDEX(J3)-1)+2)
            SYMDELTA(3*(J3-1)+3)=SYMDELTA(3*(J3-1)+3)+NEWQ(3*(SYMINDEX(J3)-1)+3)
         ENDDO
      ENDDO

C     NEWSTEPLENGTH=SUM(SYMDELTA(1:3*NATOMS)**2)
C     SYMDELTA(1:3*NATOMS)=SYMDELTA(1:3*NATOMS)*SQRT(STEPLENGTH/NEWSTEPLENGTH)
C
C  Maintain CofM distance in symmetry adapted step.
C
      SYMDELTA(1:3*NATOMS)=SYMDELTA(1:3*NATOMS)/(1+NSYMREM)

      DO J2=1,NATOMS
         ODIST1=SUM(COORDSO(3*(J2-1)+1:3*(J2-1)+3,NP)**2)
         ODIST2=SUM((COORDSO(3*(J2-1)+1:3*(J2-1)+3,NP)+SYMDELTA(3*(J2-1)+1:3*(J2-1)+3))**2)
         IF (ODIST2.NE.0.0D0) COORDS(3*(J2-1)+1:3*(J2-1)+3,NP)=(COORDSO(3*(J2-1)+1:3*(J2-1)+3,NP)
     &        +SYMDELTA(3*(J2-1)+1:3*(J2-1)+3)) * SQRT(ODIST1/ODIST2)
      ENDDO

!     DO J2=1,NATOMS
!        DO J3=J2+1,NATOMS
!           DUMMY=SQRT((COORDS(3*(J3-1)+1,NP)-COORDS(3*(J2-1)+1,NP))**2+
!    &                 (COORDS(3*(J3-1)+2,NP)-COORDS(3*(J2-1)+2,NP))**2+
!    &                 (COORDS(3*(J3-1)+3,NP)-COORDS(3*(J2-1)+3,NP))**2)
!           IF (DUMMY.LT.0.1D0) WRITE(MYUNIT,'(A,I6,A,I6,A,G20.10)') 'WARNING *** atoms ',J3,' and ',J2,
!    &             ' distance=',DUMMY
!        ENDDO
!     ENDDO

!     WRITE(789,'(I6)') NATOMS
!     WRITE(789,'(A)') 'after symmetrising step'
!     WRITE(789,'(A3,3G20.10)') ('LA ',COORDS(3*(J2-1)+1,NP),COORDS(3*(J2-1)+2,NP),
!    &                                 COORDS(3*(J2-1)+3,NP),J2=1,NATOMS-NCORE(NP))
!     WRITE(789,'(A3,3G20.10)') ('LA ',COORDS(3*(J2-1)+1,NP),COORDS(3*(J2-1)+2,NP),
!    &                                 COORDS(3*(J2-1)+3,NP),J2=NATOMS-NCORE(NP)+1,NATOMS)
!     CLOSE(789)
    
      RETURN
      END
