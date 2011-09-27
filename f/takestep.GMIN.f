C   GMIN: A program for finding global minima
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
C
      SUBROUTINE TAKESTEP(NP)
      USE commons
      IMPLICIT NONE

      DOUBLE PRECISION DPRAND, RANDOM, XMASS, YMASS, ZMASS, LOCALSTEP, DUMMY2, CDIST(NATOMS), RDOTN, XL, YL, ZL,
     1                 DIST(3*NATOMS), DMAX, VMAX, VMIN, VMAX2, EXPDECAY(NATOMS), XC, YC, ZC, ANGLE, COST, SINT,
     2                 THETA, PHI, PI, DUMMY, CMDIST(NATOMS), CMMAX, RANDOMX, RANDOMY, RANDOMZ, RVEC(3), TX, TY, TZ,
     &                 DELX, DELY, DELZ, SLENGTH, RPROJ, X(3*NATOMS), DUMMYGRAD(3*NATOMS), DUMMYE, THETA2, FRPISQ, OPOTEL
      PARAMETER (PI=3.141592654D0)
      INTEGER J1, J2, JMAX, NP, J3, JMAX2, RANATOM, INDEXD(NATOMS), NCOOPDONE, ISTAT, NDUMMY, NTRIES, NTRIESMAX

C
C  Calling CENTRE if NORESET is .TRUE. can lead to problems with COORDSO containing an atom
C  outside the permitted radius. Then it may be impossible to take a step that keeps all the
C  atoms inside.
C
      FRPISQ = 4.D0*PI*PI
      NTRIESMAX=100

C     IF (CENT.AND.(.NOT.SEEDT)) CALL CENTRE2(COORDS(1:3*NATOMS,NP)) ! COORDS might have been shifted by symmetry

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

C	csw34> For the below loop, J1 is the atom number, and J2 is set to point to 3*J1, meaning
C	the z-coordinate of that atom. The y coordinate is then J2-1 and x is J2-2. 
C
      DO J1=1,NATOMS-NSEED
         IF ((COREFRAC.EQ.0.0D0).AND.(J1.GT.NATOMS-NCORE(NP))) CYCLE ! no point taking a zero step
         NTRIES=0
10       J2=3*J1
         IF ((NCORE(NP).GT.0).AND.(J1.GT.NATOMS-NCORE(NP))) THEN
            LOCALSTEP=STEP(NP)*COREFRAC ! smaller step for core
         ELSE
            LOCALSTEP=STEP(NP)
         ENDIF
C
C  Angular move block.
C  If NORESET is .TRUE. then VAT won;t be set, so we should skip this block.
C
!        IF (J1.EQ.JMAX) WRITE(MYUNIT,'(A,I6,4F15.5)') 'JMAX,VAT,ASTEP(NP),VMIN,prod=',JMAX,VAT(J1,NP), 
!    &                                    ASTEP(NP),VMIN,ASTEP(NP)*VMIN
         ELSE IF (.NOT.FIXD) THEN
C
C  Tried changing to scale steps according to distance from CM. Maximum
C  allowed shift is linear with this distance. Worse.
C  Now try moving atoms according to how strongly bound they are.
C
           RANDOM=DPRAND()
           IF ((VMIN-VMAX.EQ.0.0D0).OR.(EFAC.EQ.0.0D0)) THEN
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
C
      ENDDO

      IF (FIXD) CALL HSMOVE(COORDS(1:3*NATOMS,1:NPAR),NP,NHSMOVE)
C
C  Preserve centre of mass if required.
C
      IF (CENT.AND.(.NOT.SEEDT)) CALL CENTRE2(COORDS(1:3*NATOMS,NP))
      IF (FIXCOM.AND.(.NOT.SEEDT)) CALL CENTRECOM(COORDS(:3*NATOMS,NP))

      RETURN
      END
C

