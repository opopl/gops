!   NEB module is an implementati on of the nudged elastic band method for performing double-ended pathway searches.
!   Copyright (C) 2003-2006 Semen A. Trygubenko and David J. Wales
!   This file is part of NEB module. NEB module is part of OPTIM.
!
!   OPTIM is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   OPTIM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
SUBROUTINE CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
USE KEY, ONLY: FROZEN, FREEZE, NREPI, NREPJ, NNREPULSIVE, DEBUG, &
  &            NCONSTRAINT, CONI, CONJ, INTCONSTRAINTDEL, CONDISTREF, INTCONSTRAINTREP, CONDISTREFLOCAL, &
  &            CONACTIVE, INTCONSTRAINREPCUT, NREPCUT,INTIMAGE
USE COMMONS, ONLY: NATOMS, NOPT
USE PORFUNCS
IMPLICIT NONE
           
INTEGER :: J1,J2,NI2,NI1,NJ2,NJ1,NMAXINT,NMININT,NREPINT(INTIMAGE+2),ISTAT
DOUBLE PRECISION :: ECON, EREP, ETOTAL, RMS
DOUBLE PRECISION R1AX,R1AY,R1AZ,R2AX,R2AY,R2AZ,R1BX,R1BY,R1BZ,R2BX,R2BY,R2BZ,D2,D1
DOUBLE PRECISION G1(3),G2(3),CONCUT,DINT,G1INT(3),G2INT(3)
DOUBLE PRECISION DUMMY, REPGRAD(3), INTCONST, D12, DSQ2, DSQ1, DSQI
DOUBLE PRECISION CONE(INTIMAGE+2), REPE(INTIMAGE+2),MAXINT,MININT,REPEINT(INTIMAGE+2),RMSIM(INTIMAGE+2)
LOGICAL NOINT
DOUBLE PRECISION XYZ(NOPT*(INTIMAGE+2)), GGG(NOPT*(INTIMAGE+2)), EEE(INTIMAGE+2)
LOGICAL IMGFREEZE(INTIMAGE)

EEE(1:INTIMAGE+2)=0.0D0
CONE(1:INTIMAGE+2)=0.0D0
REPE(1:INTIMAGE+2)=0.0D0
REPEINT(1:INTIMAGE+2)=0.0D0
NREPINT(1:INTIMAGE+2)=0
GGG(1:NOPT*(INTIMAGE+2))=0.0D0
ECON=0.0D0; EREP=0.0D0

!
!  Constraint energy and forces.
!
DO J2=1,NCONSTRAINT
   IF (.NOT.CONACTIVE(J2)) CYCLE
!
   CONCUT=0.2D0*CONDISTREF(J2)
!
! For J1 we consider the line segment between image J1-1 and J1.
! There are INTIMAGE+1 line segments in total, with an energy contribution
! and corresponding gradient terms for each. 
! A and B refer to atoms, 2 refers to image J1.
!
!  DO J1=2,INTIMAGE+1
   DO J1=1,INTIMAGE+2  ! checking for zero!
      NI1=NOPT*(J1-1)+3*(CONI(J2)-1)
      NJ1=NOPT*(J1-1)+3*(CONJ(J2)-1)
      R2AX=XYZ(NI1+1); R2AY=XYZ(NI1+2); R2AZ=XYZ(NI1+3)
      R2BX=XYZ(NJ1+1); R2BY=XYZ(NJ1+2); R2BZ=XYZ(NJ1+3)
      D2=SQRT((R2AX-R2BX)**2+(R2AY-R2BY)**2+(R2AZ-R2BZ)**2)
      IF ((ABS(D2-CONDISTREFLOCAL(J2)).GT.CONCUT)) THEN 
         DUMMY=D2-CONDISTREFLOCAL(J2)  
         G2(1)=(R2AX-R2BX)/D2
         G2(2)=(R2AY-R2BY)/D2
         G2(3)=(R2AZ-R2BZ)/D2
         REPGRAD(1:3)=2*INTCONSTRAINTDEL*((DUMMY/CONCUT)**2-1.0D0)*DUMMY*G2(1:3)
         DUMMY=INTCONSTRAINTDEL*(DUMMY**2-CONCUT**2)**2/(2.0D0*CONCUT**2)
         EEE(J1)=EEE(J1)  +DUMMY
         ECON=ECON        +DUMMY
         CONE(J1)=CONE(J1)+DUMMY
         GGG(NI1+1:NI1+3)=GGG(NI1+1:NI1+3)+REPGRAD(1:3)
         GGG(NJ1+1:NJ1+3)=GGG(NJ1+1:NJ1+3)-REPGRAD(1:3)
      ENDIF
   ENDDO
ENDDO

GGG(1:NOPT)=0.0D0                            ! can delete when loop range above changes
GGG(NOPT*(INTIMAGE+1)+1:NOPT*(INTIMAGE+2))=0.0D0 ! can delete when loop range above changes

! INTCONST=INTCONSTRAINREPCUT**13

DO J2=1,NNREPULSIVE
!  INTCONST=NREPCUT(J2)**13
   INTCONST=NREPCUT(J2)**3
!  DO J1=2,INTIMAGE+2
   DO J1=1,INTIMAGE+2 ! can change when zero energies are confirmed for end images
      NI2=NOPT*(J1-1)+3*(NREPI(J2)-1)
      NJ2=NOPT*(J1-1)+3*(NREPJ(J2)-1)
      R2AX=XYZ(NI2+1); R2AY=XYZ(NI2+2); R2AZ=XYZ(NI2+3)
      R2BX=XYZ(NJ2+1); R2BY=XYZ(NJ2+2); R2BZ=XYZ(NJ2+3)
      D2=SQRT((R2AX-R2BX)**2+(R2AY-R2BY)**2+(R2AZ-R2BZ)**2)
      IF (D2.LT.NREPCUT(J2)) THEN ! term for image J1
!        D12=D2**12
         D12=D2**2
!        DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(12.0D0*D2-13.0D0*NREPCUT(J2))/INTCONST)
         DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(2.0D0*D2-3.0D0*NREPCUT(J2))/INTCONST)
         EEE(J1)=EEE(J1)+DUMMY
         REPE(J1)=REPE(J1)+DUMMY
         EREP=EREP+DUMMY
!        DUMMY=-12.0D0*INTCONSTRAINTREP*(1.0D0/(D2*D12)-1.0D0/INTCONST)
         DUMMY=-2.0D0*INTCONSTRAINTREP*(1.0D0/(D2*D12)-1.0D0/INTCONST)
         G2(1)=(R2AX-R2BX)/D2
         G2(2)=(R2AY-R2BY)/D2
         G2(3)=(R2AZ-R2BZ)/D2
         REPGRAD(1:3)=DUMMY*G2(1:3)
         GGG(NI2+1:NI2+3)=GGG(NI2+1:NI2+3)+REPGRAD(1:3)
         GGG(NJ2+1:NJ2+3)=GGG(NJ2+1:NJ2+3)-REPGRAD(1:3)
      ENDIF
!
! For internal minima we are counting edges. 
! Edge J1 is between images J1-1 and J1, starting from J1=2.
! Energy contributions are shared evenly, except for
! edge 1, which is assigned to image 2, and edge INTIMAGE+1, which
! is assigned to image INTIMAGE+1. Gradients are set to zero for
! the end images.
!
      IF (J1.EQ.1) CYCLE
      NI1=NOPT*(J1-2)+3*(NREPI(J2)-1)
      NJ1=NOPT*(J1-2)+3*(NREPJ(J2)-1)
      R1AX=XYZ(NI1+1); R1AY=XYZ(NI1+2); R1AZ=XYZ(NI1+3)
      R1BX=XYZ(NJ1+1); R1BY=XYZ(NJ1+2); R1BZ=XYZ(NJ1+3)
      CALL MINMAXD2R(R1AX,R1AY,R1AZ,R2AX,R2AY,R2AZ,R1BX,R1BY,R1BZ,R2BX,R2BY,R2BZ, &
  &                 D2,D1,DINT,DSQ2,DSQ1,DSQI,G1,G2,G1INT,G2INT,NOINT,.FALSE.,NREPCUT(J2))
      IF ((.NOT.NOINT).AND.(DINT.LT.NREPCUT(J2))) THEN
!        D12=DSQI**6
         D12=DSQI
!        DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(12.0D0*DINT-13.0D0*NREPCUT(J2))/INTCONST)
         DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(2.0D0*DINT-3.0D0*NREPCUT(J2))/INTCONST)
         IF (J1.EQ.2) THEN
            EEE(J1)=EEE(J1)+DUMMY
            REPEINT(J1)=REPEINT(J1)+DUMMY
            NREPINT(J1)=NREPINT(J1)+1
         ELSE IF (J1.LT.INTIMAGE+2) THEN
            EEE(J1)=EEE(J1)+DUMMY/2.0D0
            EEE(J1-1)=EEE(J1-1)+DUMMY/2.0D0
            REPEINT(J1)=REPEINT(J1)+DUMMY/2.0D0
            REPEINT(J1-1)=REPEINT(J1-1)+DUMMY/2.0D0
            NREPINT(J1)=NREPINT(J1)+1
            NREPINT(J1-1)=NREPINT(J1-1)+1
         ELSE IF (J1.EQ.INTIMAGE+2) THEN
            EEE(J1-1)=EEE(J1-1)+DUMMY
            REPEINT(J1-1)=REPEINT(J1-1)+DUMMY
            NREPINT(J1-1)=NREPINT(J1-1)+1
         ENDIF
         EREP=EREP+DUMMY
!        DUMMY=-12.0D0*INTCONSTRAINTREP*(1.0D0/(DINT*D12)-1.0D0/INTCONST)
         DUMMY=-2.0D0*INTCONSTRAINTREP*(1.0D0/(DINT*D12)-1.0D0/INTCONST)
         REPGRAD(1:3)=DUMMY*G1INT(1:3)
!        PRINT '(A,4I6,2G15.5)','in1 J1,J2,REPI,REPJ,REPGRAD,NREPCUT=',J1,J2,NREPI(J2),NREPJ(J2), &
! &                              SQRT(REPGRAD(1)**2+REPGRAD(2)**2+REPGRAD(3)**2),NREPCUT(J2)
!
! Gradient contributions for image J1-1
!
         GGG(NI1+1:NI1+3)=GGG(NI1+1:NI1+3)+REPGRAD(1:3)
         GGG(NJ1+1:NJ1+3)=GGG(NJ1+1:NJ1+3)-REPGRAD(1:3)
         REPGRAD(1:3)=DUMMY*G2INT(1:3)
!        PRINT '(A,4I6,2G15.5)','in1 J1,J2,REPI,REPJ,REPGRAD,NREPCUT=',J1,J2,NREPI(J2),NREPJ(J2), &
! &                              SQRT(REPGRAD(1)**2+REPGRAD(2)**2+REPGRAD(3)**2),NREPCUT(J2)
!
! Gradient contributions for image J1
!
         GGG(NI2+1:NI2+3)=GGG(NI2+1:NI2+3)+REPGRAD(1:3)
         GGG(NJ2+1:NJ2+3)=GGG(NJ2+1:NJ2+3)-REPGRAD(1:3)
      ENDIF
   ENDDO
ENDDO
!
! Set gradients on frozen atoms to zero.
!
IF (FREEZE) THEN
   DO J1=2,INTIMAGE+1  
      DO J2=1,NATOMS
         IF (FROZEN(J2)) THEN
            GGG(NOPT*(J1-1)+3*(J2-1)+1)=0.0D0
            GGG(NOPT*(J1-1)+3*(J2-1)+2)=0.0D0
            GGG(NOPT*(J1-1)+3*(J2-1)+3)=0.0D0
         ENDIF
      ENDDO
   ENDDO
ENDIF
!
! Set gradients to zero for start and finish images.
!
GGG(1:NOPT)=0.0D0
GGG((INTIMAGE+1)*NOPT+1:(INTIMAGE+2)*NOPT)=0.0D0
RMS=0.0D0
DO J1=2,INTIMAGE+1
   RMSIM(J1)=0.0D0
   DO J2=1,NOPT
      RMS=RMS+GGG(NOPT*(J1-1)+J2)**2
      RMSIM(J1)=RMSIM(J1)+GGG(NOPT*(J1-1)+J2)**2
   ENDDO
   RMSIM(J1)=SQRT(RMSIM(J1)/NOPT)
ENDDO
IF (INTIMAGE.NE.0) THEN
   RMS=SQRT(RMS/(NOPT*INTIMAGE))
ENDIF
!
! For INTIMAGE images there are INTIMAGE+2 replicas including the end points,
! and INTIMAGE+1 line segements, with associated energies stored in EEE(2:INTIMAGE+2)
!
ETOTAL=0.0D0
MAXINT=-1.0D100
MININT=1.0D100
DO J1=2,INTIMAGE+1
   ETOTAL=ETOTAL+EEE(J1)
   PRINT '(A,I6,A,3G20.10)',' congrad> con/rep/RMS image ',J1,' ',CONE(J1),REPE(J1),RMSIM(J1)
   IF (REPEINT(J1).LT.MININT) THEN
      MININT=REPEINT(J1)
      NMININT=J1
   ENDIF
   IF (REPE(J1).GT.MAXINT) THEN
      MAXINT=REPE(J1)
      NMAXINT=J1
   ENDIF
ENDDO
IF (DEBUG) PRINT '(A,G20.10,A,2I6)',' congrad> largest  internal energy=',MAXINT,' for image ',NMAXINT
IF (DEBUG) PRINT '(A,G20.10,A,2I6)',' congrad> smallest internal energy=',MININT,' for image ',NMININT

END SUBROUTINE CONGRAD

SUBROUTINE MINMAXD2(R1AX,R1AY,R1AZ,R2AX,R2AY,R2AZ,R1BX,R1BY,R1BZ,R2BX,R2BY,R2BZ, &
  &                 D2,D1,DINT,G1,G2,G1INT,G2INT,NOINT,DEBUG)
IMPLICIT NONE
DOUBLE PRECISION R1AX,R1AY,R1AZ,R2AX,R2AY,R2AZ,R1BX,R1BY,R1BZ,R2BX,R2BY,R2BZ,D2,D1,DINT
DOUBLE PRECISION G1(3),G2(3),G1INT(3),G2INT(3)
DOUBLE PRECISION DSQ2, DSQ1, DSQI, r1apr2bmr2amr1bsq, r1amr1bsq, r2amr2bsq
DOUBLE PRECISION r1amr1bdr2amr2b, r1amr1bdr2amr2bsq, DUMMY
LOGICAL NOINT, DEBUG
!
! Squared distance between atoms A and B for theta=0 - distance in image 2
!
DSQ2=r2ax**2 + r2ay**2 + r2az**2 + r2bx**2 + r2by**2 + r2bz**2 - 2*(r2ax*r2bx + r2ay*r2by + r2az*r2bz)
!
! Squared distance between atoms A and B for theta=Pi/2 - distance in image 1
!
DSQ1=r1ax**2 + r1ay**2 + r1az**2 + r1bx**2 + r1by**2 + r1bz**2 - 2*(r1ax*r1bx + r1ay*r1by + r1az*r1bz)
! PRINT '(A,6F15.10)','R1AX,R1AY,R1AZ,R2AX,R2AY,R2AZ=',R1AX,R1AY,R1AZ,R2AX,R2AY,R2AZ
! PRINT '(A,6F15.10)','R1BX,R1BY,R1BZ,R2BX,R2BY,R2BZ=',R1BX,R1BY,R1BZ,R2BX,R2BY,R2BZ
!
! Is there an internal extremum?
!
r1apr2bmr2amr1bsq=(r1ax-r1bx-r2ax+r2bx)**2+(r1ay-r1by-r2ay+r2by)**2+(r1az-r1bz-r2az+r2bz)**2
IF (r1apr2bmr2amr1bsq.EQ.0.0D0) THEN
   DUMMY=2.0D0 ! just to skip the internal extremum part
ELSE
   DUMMY=((r1ax-r1bx)*(r1ax-r1bx-r2ax+r2bx)+(r1ay-r1by)*(r1ay-r1by-r2ay+r2by)+(r1az-r1bz)*(r1az-r1bz-r2az+r2bz))/r1apr2bmr2amr1bsq
ENDIF
NOINT=.TRUE.
IF ((DUMMY.GT.0.0D0).AND.(DUMMY.LT.1.0D0)) NOINT=.FALSE.
G2(1:3)=0.0D0
G1(1:3)=0.0D0
G1INT(1:3)=0.0D0
G2INT(1:3)=0.0D0
D2=SQRT(DSQ2)
D1=SQRT(DSQ1)
G2(1)=r2ax - r2bx
G2(2)=r2ay - r2by
G2(3)=r2az - r2bz
G1(1)=r1ax - r1bx
G1(2)=r1ay - r1by
G1(3)=r1az - r1bz
DSQI=1.0D10
DINT=1.0D10
IF (.NOT.NOINT) THEN
   r1amr1bdr2amr2b=(r1ax-r1bx)*(r2ax-r2bx)+(r1ay-r1by)*(r2ay-r2by)+(r1az-r1bz)*(r2az-r2bz)
   r1amr1bdr2amr2bsq=r1amr1bdr2amr2b**2
   r1amr1bsq=(r1ax - r1bx)**2 + (r1ay - r1by)**2 + (r1az - r1bz)**2
   r2amr2bsq=(r2ax - r2bx)**2 + (r2ay - r2by)**2 + (r2az - r2bz)**2
   DSQI=(-r1amr1bdr2amr2bsq + r1amr1bsq*r2amr2bsq)/r1apr2bmr2amr1bsq
   DUMMY=r1apr2bmr2amr1bsq**2
   DINT=SQRT(DSQI)
   IF (DINT.LE.0.0D0) THEN
      NOINT=.TRUE.
   ELSE
      G1INT(1)= (((r1amr1bdr2amr2bsq - r1amr1bsq*r2amr2bsq)*(r1ax - r1bx - r2ax + r2bx) + &
 &    r1apr2bmr2amr1bsq*((r1ax - r1bx)*r2amr2bsq + r1amr1bdr2amr2b*(-r2ax + r2bx))))/DUMMY
      G1INT(2)= (((r1amr1bdr2amr2bsq - r1amr1bsq*r2amr2bsq)*(r1ay - r1by - r2ay + r2by) + &
 &    r1apr2bmr2amr1bsq*((r1ay - r1by)*r2amr2bsq + r1amr1bdr2amr2b*(-r2ay + r2by))))/DUMMY
      G1INT(3)= (((r1amr1bdr2amr2bsq - r1amr1bsq*r2amr2bsq)*(r1az - r1bz - r2az + r2bz) + &
 &    r1apr2bmr2amr1bsq*((r1az - r1bz)*r2amr2bsq + r1amr1bdr2amr2b*(-r2az + r2bz))))/DUMMY

      G2INT(1)= (((r1amr1bdr2amr2bsq - r1amr1bsq*r2amr2bsq)*(r2ax - r2bx - r1ax + r1bx) + &
 &    r1apr2bmr2amr1bsq*((r2ax - r2bx)*r1amr1bsq + r1amr1bdr2amr2b*(-r1ax + r1bx))))/DUMMY
      G2INT(2)= (((r1amr1bdr2amr2bsq - r1amr1bsq*r2amr2bsq)*(r2ay - r2by - r1ay + r1by) + &
 &    r1apr2bmr2amr1bsq*((r2ay - r2by)*r1amr1bsq + r1amr1bdr2amr2b*(-r1ay + r1by))))/DUMMY
      G2INT(3)= (((r1amr1bdr2amr2bsq - r1amr1bsq*r2amr2bsq)*(r2az - r2bz - r1az + r1bz) + &
 &    r1apr2bmr2amr1bsq*((r2az - r2bz)*r1amr1bsq + r1amr1bdr2amr2b*(-r1az + r1bz))))/DUMMY
   ENDIF
ENDIF
!
! Convert derivatives of distance^2 to derivative of distance.
! We have cancelled a factor of two above and below!
!
G2(1:3)=G2(1:3)/D2
G1(1:3)=G1(1:3)/D1
IF (.NOT.NOINT) THEN
   G1INT(1:3)=G1INT(1:3)/DINT
   G2INT(1:3)=G2INT(1:3)/DINT
ENDIF

! PRINT '(A,3G12.5,L5)','D2,D1,DINT,NOINT=',D2,D1,DINT,NOINT

END SUBROUTINE MINMAXD2

SUBROUTINE MINMAXD2R(R1AX,R1AY,R1AZ,R2AX,R2AY,R2AZ,R1BX,R1BY,R1BZ,R2BX,R2BY,R2BZ, &
  &                 D2,D1,DINT,DSQ2,DSQ1,DSQI,G1,G2,G1INT,G2INT,NOINT,DEBUG,INTCONSTRAINREPCUT)
IMPLICIT NONE
DOUBLE PRECISION R1AX,R1AY,R1AZ,R2AX,R2AY,R2AZ,R1BX,R1BY,R1BZ,R2BX,R2BY,R2BZ,D2,D1,DINT
DOUBLE PRECISION G1(3),G2(3),G1INT(3),G2INT(3),INTCONSTRAINREPCUT
DOUBLE PRECISION DSQ2, DSQ1, DSQI, r1apr2bmr2amr1bsq, r1amr1bsq, r2amr2bsq
DOUBLE PRECISION r1amr1bdr2amr2b, r1amr1bdr2amr2bsq, DUMMY
LOGICAL NOINT, DEBUG
!
! Squared distance between atoms A and B for theta=0 - distance in image 2
!
DSQ2=r2ax**2 + r2ay**2 + r2az**2 + r2bx**2 + r2by**2 + r2bz**2 - 2*(r2ax*r2bx + r2ay*r2by + r2az*r2bz)
!
! Squared distance between atoms A and B for theta=Pi/2 - distance in image 1
!
DSQ1=r1ax**2 + r1ay**2 + r1az**2 + r1bx**2 + r1by**2 + r1bz**2 - 2*(r1ax*r1bx + r1ay*r1by + r1az*r1bz)
!
! Is there an internal extremum?
!
r1apr2bmr2amr1bsq=(r1ax-r1bx-r2ax+r2bx)**2+(r1ay-r1by-r2ay+r2by)**2+(r1az-r1bz-r2az+r2bz)**2
IF (r1apr2bmr2amr1bsq.EQ.0.0D0) THEN
   DUMMY=2.0D0 ! just to skip the internal solution
ELSE
   DUMMY=((r1ax-r1bx)*(r1ax-r1bx-r2ax+r2bx)+(r1ay-r1by)*(r1ay-r1by-r2ay+r2by)+(r1az-r1bz)*(r1az-r1bz-r2az+r2bz))/r1apr2bmr2amr1bsq
ENDIF
NOINT=.TRUE.
IF ((DUMMY.GT.0.0D0).AND.(DUMMY.LT.1.0D0)) NOINT=.FALSE.
G2(1:3)=0.0D0
G1(1:3)=0.0D0
G1INT(1:3)=0.0D0
G2INT(1:3)=0.0D0
D2=SQRT(DSQ2)
D1=SQRT(DSQ1)
G2(1)=r2ax - r2bx
G2(2)=r2ay - r2by
G2(3)=r2az - r2bz
G1(1)=r1ax - r1bx
G1(2)=r1ay - r1by
G1(3)=r1az - r1bz
DSQI=1.0D10
DINT=1.0D10
IF (.NOT.NOINT) THEN
   r1amr1bdr2amr2b=(r1ax-r1bx)*(r2ax-r2bx)+(r1ay-r1by)*(r2ay-r2by)+(r1az-r1bz)*(r2az-r2bz)
   r1amr1bdr2amr2bsq=r1amr1bdr2amr2b**2
   r1amr1bsq=(r1ax - r1bx)**2 + (r1ay - r1by)**2 + (r1az - r1bz)**2
   r2amr2bsq=(r2ax - r2bx)**2 + (r2ay - r2by)**2 + (r2az - r2bz)**2
   DSQI=MAX((-r1amr1bdr2amr2bsq + r1amr1bsq*r2amr2bsq)/r1apr2bmr2amr1bsq,0.0D0)
   DUMMY=r1apr2bmr2amr1bsq**2
   DINT=SQRT(DSQI)
   IF (DINT.LE.0.0D0) THEN
      NOINT=.TRUE.
   ELSEIF (DINT.LE.INTCONSTRAINREPCUT) THEN ! skip otherwise
      G1INT(1)= (((r1amr1bdr2amr2bsq - r1amr1bsq*r2amr2bsq)*(r1ax - r1bx - r2ax + r2bx) + &
 &    r1apr2bmr2amr1bsq*((r1ax - r1bx)*r2amr2bsq + r1amr1bdr2amr2b*(-r2ax + r2bx))))/DUMMY
      G1INT(2)= (((r1amr1bdr2amr2bsq - r1amr1bsq*r2amr2bsq)*(r1ay - r1by - r2ay + r2by) + &
 &    r1apr2bmr2amr1bsq*((r1ay - r1by)*r2amr2bsq + r1amr1bdr2amr2b*(-r2ay + r2by))))/DUMMY
      G1INT(3)= (((r1amr1bdr2amr2bsq - r1amr1bsq*r2amr2bsq)*(r1az - r1bz - r2az + r2bz) + &
 &    r1apr2bmr2amr1bsq*((r1az - r1bz)*r2amr2bsq + r1amr1bdr2amr2b*(-r2az + r2bz))))/DUMMY

      G2INT(1)= (((r1amr1bdr2amr2bsq - r1amr1bsq*r2amr2bsq)*(r2ax - r2bx - r1ax + r1bx) + &
 &    r1apr2bmr2amr1bsq*((r2ax - r2bx)*r1amr1bsq + r1amr1bdr2amr2b*(-r1ax + r1bx))))/DUMMY
      G2INT(2)= (((r1amr1bdr2amr2bsq - r1amr1bsq*r2amr2bsq)*(r2ay - r2by - r1ay + r1by) + &
 &    r1apr2bmr2amr1bsq*((r2ay - r2by)*r1amr1bsq + r1amr1bdr2amr2b*(-r1ay + r1by))))/DUMMY
      G2INT(3)= (((r1amr1bdr2amr2bsq - r1amr1bsq*r2amr2bsq)*(r2az - r2bz - r1az + r1bz) + &
 &    r1apr2bmr2amr1bsq*((r2az - r2bz)*r1amr1bsq + r1amr1bdr2amr2b*(-r1az + r1bz))))/DUMMY
   ENDIF
ENDIF
!
! Convert derivatives of distance^2 to derivative of distance.
! We have cancelled a factor of two above and below!
!
G2(1:3)=G2(1:3)/D2
G1(1:3)=G1(1:3)/D1
IF (.NOT.NOINT) THEN
!  IF (DINT.EQ.0.0D0) THEN
!     PRINT '(A,G20.10)','minmaxd2r> ERROR *** DINT=',DINT
!     PRINT *,'original dummy=',((r1ax-r1bx)*(r1ax-r1bx-r2ax+r2bx)+ &
! &        (r1ay-r1by)*(r1ay-r1by-r2ay+r2by)+(r1az-r1bz)*(r1az-r1bz-r2az+r2bz))/r1apr2bmr2amr1bsq
!     PRINT *,'r1amr1bdr2amr2b=',r1amr1bdr2amr2b
!     PRINT *,'r1amr1bdr2amr2bsq=',r1amr1bdr2amr2bsq
!     PRINT *,'r1amr1bsq=',r1amr1bsq
!     PRINT *,'r2amr2bsq=',r2amr2bsq
!     PRINT *,'DSQI=',DSQI
!     PRINT *,'DUMMY=',DUMMY
!     PRINT *,'G1INT=',G1INT(1:3)
!     PRINT *,'G2INT=',G2INT(1:3)
!     PRINT *,'R1AX,R1AY,R1AZ=',R1AX,R1AY,R1AZ
!     PRINT *,'R2AX,R2AY,R2AZ=',R2AX,R2AY,R2AZ
!     PRINT *,'R1BX,R1BY,R1BZ=',R1BX,R1BY,R1BZ
!     PRINT *,'R2BX,R2BY,R2BZ=',R2BX,R2BY,R2BZ
!     STOP
!  ENDIF
   G1INT(1:3)=G1INT(1:3)/DINT
   G2INT(1:3)=G2INT(1:3)/DINT
ENDIF

END SUBROUTINE MINMAXD2R

!
! This version of congrad tests for an internal minimum in the
! constraint distances as well as the repulsions.
!
SUBROUTINE CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
USE KEY, ONLY: FROZEN, FREEZE, NREPI, NREPJ, NNREPULSIVE, DEBUG, &
  &            NCONSTRAINT, CONI, CONJ, INTCONSTRAINTDEL, CONDISTREF, INTCONSTRAINTREP, CONDISTREFLOCAL, &
  &            CONACTIVE, INTCONSTRAINREPCUT, NREPCUT,FREEZENODEST, INTIMAGE
USE COMMONS, ONLY: NATOMS, NOPT
IMPLICIT NONE
           
INTEGER :: J1,J2,NI2,NI1,NJ2,NJ1,NMAXINT,NMININT,NCONINT(INTIMAGE+2),NREPINT(INTIMAGE+2)
DOUBLE PRECISION :: ECON, EREP, ETOTAL, RMS
DOUBLE PRECISION R1AX,R1AY,R1AZ,R2AX,R2AY,R2AZ,R1BX,R1BY,R1BZ,R2BX,R2BY,R2BZ,D2,D1
DOUBLE PRECISION G1(3),G2(3),CONCUT,DINT,G1INT(3),G2INT(3)
DOUBLE PRECISION DUMMY, REPGRAD(3), INTCONST, D12, DSQ2, DSQ1, DSQI
DOUBLE PRECISION CONE(INTIMAGE+2), REPE(INTIMAGE+2),MAXINT,MININT,REPEINT(INTIMAGE+2),CONEINT(INTIMAGE+2),RMSIMAGE(INTIMAGE+2)
LOGICAL NOINT, LPRINT
DOUBLE PRECISION XYZ(NOPT*(INTIMAGE+2)), GGG(NOPT*(INTIMAGE+2)), EEE(INTIMAGE+2)
LOGICAL IMGFREEZE(INTIMAGE)

EEE(1:INTIMAGE+2)=0.0D0
CONE(1:INTIMAGE+2)=0.0D0
REPE(1:INTIMAGE+2)=0.0D0
NCONINT(1:INTIMAGE+2)=0
NREPINT(1:INTIMAGE+2)=0
REPEINT(1:INTIMAGE+2)=0.0D0
CONEINT(1:INTIMAGE+2)=0.0D0
GGG(1:NOPT*(INTIMAGE+2))=0.0D0
ECON=0.0D0; EREP=0.0D0
LPRINT=.TRUE.
LPRINT=.FALSE.
!
!  Constraint energy and forces.
!
! For J1 we consider the line segment between image J1-1 and J1.
! There are INTIMAGE+1 line segments in total, with an energy contribution
! and corresponding gradient terms for each. 
! A and B refer to atoms, 1 and 2 to images J1-1 and J1 corresponding to J1-2 and J1-1 below.
!
! IMGFREEZE(1:INTIMAGE) refers to the images excluding end points!
!
DO J2=1,NCONSTRAINT
   IF (.NOT.CONACTIVE(J2)) CYCLE
   DO J1=2,INTIMAGE+2
      IF (FREEZENODEST) THEN ! IMGFREEZE is not allocated otherwise!
         IF (J1.EQ.2) THEN
            IF (IMGFREEZE(1)) CYCLE
         ELSE IF (J1.EQ.INTIMAGE+2) THEN
            IF (IMGFREEZE(INTIMAGE)) CYCLE
         ELSE
            IF (IMGFREEZE(J1-2).AND.IMGFREEZE(J1-1)) CYCLE
         ENDIF
      ENDIF
      NI1=NOPT*(J1-2)+3*(CONI(J2)-1)
      NI2=NOPT*(J1-1)+3*(CONI(J2)-1)
      NJ1=NOPT*(J1-2)+3*(CONJ(J2)-1)
      NJ2=NOPT*(J1-1)+3*(CONJ(J2)-1)
      R1AX=XYZ(NI1+1); R1AY=XYZ(NI1+2); R1AZ=XYZ(NI1+3)
      R1BX=XYZ(NJ1+1); R1BY=XYZ(NJ1+2); R1BZ=XYZ(NJ1+3)
      R2AX=XYZ(NI2+1); R2AY=XYZ(NI2+2); R2AZ=XYZ(NI2+3)
      R2BX=XYZ(NJ2+1); R2BY=XYZ(NJ2+2); R2BZ=XYZ(NJ2+3)
      CALL MINMAXD2(R1AX,R1AY,R1AZ,R2AX,R2AY,R2AZ,R1BX,R1BY,R1BZ,R2BX,R2BY,R2BZ, &
  &                 D2,D1,DINT,G1,G2,G1INT,G2INT,NOINT,.FALSE.)
!
! Need to include both D2 and D1 contributions if they are both outside tolerance.
! Otherwise we get discontinuities if they are very close and swap over.
!
      CONCUT=CONDISTREF(J2)*0.2D0
!
! terms for image J1 - non-zero derivatives only for J1. D2 is the distance for image J1.
!
!     IF (LPRINT) PRINT '(A,I6,5G15.5)', &
! &       'J1,D2,D1,DINT,MIN diff,CONCUT=',J1,D2,D1,DINT,ABS(D2-CONDISTREFLOCAL(J2)),CONCUT
      IF ((ABS(D2-CONDISTREFLOCAL(J2)).GT.CONCUT).AND.(J1.LT.INTIMAGE+2)) THEN 
         DUMMY=D2-CONDISTREFLOCAL(J2)  
         REPGRAD(1:3)=2*INTCONSTRAINTDEL*((DUMMY/CONCUT)**2-1.0D0)*DUMMY*G2(1:3)
         DUMMY=INTCONSTRAINTDEL*(DUMMY**2-CONCUT**2)**2/(2.0D0*CONCUT**2)
         EEE(J1)=EEE(J1)+DUMMY
         CONE(J1)=CONE(J1)+DUMMY
         ECON=ECON      +DUMMY
!        IF (LPRINT) PRINT '(A,4I6,G15.5)','min J1,J2,CONI,CONJ,REPGRAD=',J1,J2,CONI(J2),CONJ(J2), &
! &                              SQRT(REPGRAD(1)**2+REPGRAD(2)**2+REPGRAD(3)**2)
         GGG(NI2+1:NI2+3)=GGG(NI2+1:NI2+3)+REPGRAD(1:3)
         GGG(NJ2+1:NJ2+3)=GGG(NJ2+1:NJ2+3)-REPGRAD(1:3)
      ENDIF
!
! Don't add energy contributions to EEE(2) from D1, since the gradients are non-zero only for image 1.
!
! terms for image J1-1 - non-zero derivatives only for J1-1. D1 is the distance for image J1-1.
!
!     IF (LPRINT) PRINT '(A,I6,5G15.5)', &
! &       'J1,D2,D1,DINT,MAX diff,CONCUT=',J1,D2,D1,DINT,ABS(D1-CONDISTREFLOCAL(J2)),CONCUT
      IF ((ABS(D1-CONDISTREFLOCAL(J2)).GT.CONCUT).AND.(J1.GT.2)) THEN  
         DUMMY=D1-CONDISTREFLOCAL(J2)  
         REPGRAD(1:3)=2*INTCONSTRAINTDEL*((DUMMY/CONCUT)**2-1.0D0)*DUMMY*G1(1:3)
         DUMMY=INTCONSTRAINTDEL*(DUMMY**2-CONCUT**2)**2/(2.0D0*CONCUT**2)
         EEE(J1-1)=EEE(J1-1)+DUMMY
         CONE(J1-1)=CONE(J1-1)+DUMMY
         ECON=ECON      +DUMMY
!        IF (LPRINT) PRINT '(A,4I6,G15.5)','max J1,J2,CONI,CONJ,REPGRAD=',J1,J2,CONI(J2),CONJ(J2), &
! &         SQRT(REPGRAD(1)**2+REPGRAD(2)**2+REPGRAD(3)**2)
         GGG(NI1+1:NI1+3)=GGG(NI1+1:NI1+3)+REPGRAD(1:3)
         GGG(NJ1+1:NJ1+3)=GGG(NJ1+1:NJ1+3)-REPGRAD(1:3)
      ENDIF
      IF ((.NOT.NOINT).AND.(ABS(DINT-CONDISTREFLOCAL(J2)).GT.CONCUT)) THEN
         DUMMY=DINT-CONDISTREFLOCAL(J2)  
         REPGRAD(1:3)=2*INTCONSTRAINTDEL*((DUMMY/CONCUT)**2-1.0D0)*DUMMY*G1INT(1:3)
         GGG(NI1+1:NI1+3)=GGG(NI1+1:NI1+3)+REPGRAD(1:3)
         GGG(NJ1+1:NJ1+3)=GGG(NJ1+1:NJ1+3)-REPGRAD(1:3)
         REPGRAD(1:3)=2*INTCONSTRAINTDEL*((DUMMY/CONCUT)**2-1.0D0)*DUMMY*G2INT(1:3)
         DUMMY=INTCONSTRAINTDEL*(DUMMY**2-CONCUT**2)**2/(2.0D0*CONCUT**2)
         ECON=ECON+DUMMY
         IF (J1.EQ.2) THEN
            EEE(J1)=EEE(J1)+DUMMY
            CONEINT(J1)=CONEINT(J1)+DUMMY
            NCONINT(J1)=NCONINT(J1)+1
         ELSE IF (J1.LT.INTIMAGE+2) THEN
            EEE(J1)=EEE(J1)+DUMMY/2.0D0
            EEE(J1-1)=EEE(J1-1)+DUMMY/2.0D0
            CONEINT(J1)=CONEINT(J1)+DUMMY/2.0D0
            CONEINT(J1-1)=CONEINT(J1-1)+DUMMY/2.0D0
            NCONINT(J1)=NCONINT(J1)+1
            NCONINT(J1-1)=NCONINT(J1-1)+1
         ELSE IF (J1.EQ.INTIMAGE+2) THEN
            EEE(J1-1)=EEE(J1-1)+DUMMY
            CONEINT(J1-1)=CONEINT(J1-1)+DUMMY
            NCONINT(J1-1)=NCONINT(J1-1)+1
         ENDIF
!        PRINT '(A,4I6,G15.5)','in2 J1,J2,CONI,CONJ,REPGRAD=',J1,J2,CONI(J2),CONJ(J2), &
! &                              SQRT(REPGRAD(1)**2+REPGRAD(2)**2+REPGRAD(3)**2)
         GGG(NI2+1:NI2+3)=GGG(NI2+1:NI2+3)+REPGRAD(1:3)
         GGG(NJ2+1:NJ2+3)=GGG(NJ2+1:NJ2+3)-REPGRAD(1:3)
      ENDIF
   ENDDO
ENDDO

! INTCONST=INTCONSTRAINREPCUT**13

DO J2=1,NNREPULSIVE
!  INTCONST=NREPCUT(J2)**13
   INTCONST=NREPCUT(J2)**3
   DO J1=2,INTIMAGE+2
      If (FREEZENODEST) THEN
         IF (J1.EQ.2) THEN
            IF (IMGFREEZE(1)) CYCLE
         ELSE IF (J1.EQ.INTIMAGE+2) THEN
            IF (IMGFREEZE(INTIMAGE)) CYCLE
         ELSE
            IF (IMGFREEZE(J1-2).AND.IMGFREEZE(J1-1)) CYCLE
         ENDIF
      ENDIF
      NI1=NOPT*(J1-2)+3*(NREPI(J2)-1)
      NI2=NOPT*(J1-1)+3*(NREPI(J2)-1)
      NJ1=NOPT*(J1-2)+3*(NREPJ(J2)-1)
      NJ2=NOPT*(J1-1)+3*(NREPJ(J2)-1)
      R1AX=XYZ(NI1+1); R1AY=XYZ(NI1+2); R1AZ=XYZ(NI1+3)
      R1BX=XYZ(NJ1+1); R1BY=XYZ(NJ1+2); R1BZ=XYZ(NJ1+3)
      R2AX=XYZ(NI2+1); R2AY=XYZ(NI2+2); R2AZ=XYZ(NI2+3)
      R2BX=XYZ(NJ2+1); R2BY=XYZ(NJ2+2); R2BZ=XYZ(NJ2+3)
      CALL MINMAXD2R(R1AX,R1AY,R1AZ,R2AX,R2AY,R2AZ,R1BX,R1BY,R1BZ,R2BX,R2BY,R2BZ, &
  &                 D2,D1,DINT,DSQ2,DSQ1,DSQI,G1,G2,G1INT,G2INT,NOINT,.FALSE.,NREPCUT(J2))
!     IF ((NREPI(J2).EQ.135).AND.(NREPJ(J2).EQ.192)) THEN
!        PRINT '(A,3G20.10)',' congrad2> R1AX,R1AY,R1AZ=',R1AX,R1AY,R1AZ
!        PRINT '(A,3G20.10)',' congrad2> R1BX,R1BY,R1BZ=',R1BX,R1BY,R1BZ
!        PRINT '(A,3G20.10)',' congrad2> R2AX,R2AY,R2AZ=',R2AX,R2AY,R2AZ
!        PRINT '(A,3G20.10)',' congrad2> R2BX,R2BY,R2BZ=',R2BX,R2BY,R2BZ
!        PRINT '(A,I6,A,2I6)',' congrad2> J1=',J1,' edge between images: ',J1-1,J1
!        PRINT '(A,L5,3G20.10)',' congrad2> NOINT,D2,D1,DINT=',NOINT,D2,D1,DINT
!     ENDIF
      DUMMY=0.0D0 
!
! Skip image INTIMAGE+2 - no non-zero gradients on other images and no energy contributions.
!
!     IF ((D2.LT.INTCONSTRAINREPCUT).AND.(J1.LT.INTIMAGE+2)) THEN ! terms for image J1 - non-zero derivatives only for J1
      IF ((D2.LT.NREPCUT(J2)).AND.(J1.LT.INTIMAGE+2)) THEN ! terms for image J1 - non-zero derivatives only for J1
!        D12=DSQ2**6
         D12=DSQ2
!        DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(12.0D0*D2-13.0D0*INTCONSTRAINREPCUT)/INTCONST)
!        DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(12.0D0*D2-13.0D0*NREPCUT(J2))/INTCONST)
         DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(2.0D0*D2-3.0D0*NREPCUT(J2))/INTCONST)
         EEE(J1)=EEE(J1)+DUMMY
         REPE(J1)=REPE(J1)+DUMMY
         EREP=EREP+DUMMY
!        DUMMY=-12.0D0*INTCONSTRAINTREP*(1.0D0/(D2*D12)-1.0D0/INTCONST)
         DUMMY=-2.0D0*INTCONSTRAINTREP*(1.0D0/(D2*D12)-1.0D0/INTCONST)
         REPGRAD(1:3)=DUMMY*G2(1:3)
!        PRINT '(A,4I6,G15.5)','min J1,J2,REPI,REPJ,REPGRAD=',J1,J2,NREPI(J2),NREPJ(J2), &
! &                              SQRT(REPGRAD(1)**2+REPGRAD(2)**2+REPGRAD(3)**2)
         GGG(NI2+1:NI2+3)=GGG(NI2+1:NI2+3)+REPGRAD(1:3)
         GGG(NJ2+1:NJ2+3)=GGG(NJ2+1:NJ2+3)-REPGRAD(1:3)
      ENDIF
      DUMMY=0.0D0
!
! Don't add energy contributions to EEE(2) from D1, since the gradients are non-zero only for image 1.
!
!     IF ((D1.LT.INTCONSTRAINREPCUT).AND.(J1.GT.2)) THEN ! terms for image J1-1 - non-zero derivatives only for J1-1
      IF ((D1.LT.NREPCUT(J2)).AND.(J1.GT.2)) THEN ! terms for image J1-1 - non-zero derivatives only for J1-1
!        D12=DSQ1**6
         D12=DSQ1
!        DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(12.0D0*D1-13.0D0*INTCONSTRAINREPCUT)/INTCONST)
!        DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(12.0D0*D1-13.0D0*NREPCUT(J2))/INTCONST)
         DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(2.0D0*D1-3.0D0*NREPCUT(J2))/INTCONST)
         EEE(J1-1)=EEE(J1-1)+DUMMY
         REPE(J1-1)=REPE(J1-1)+DUMMY
         EREP=EREP+DUMMY
!        DUMMY=-12.0D0*INTCONSTRAINTREP*(1.0D0/(D1*D12)-1.0D0/INTCONST)
         DUMMY=-2.0D0*INTCONSTRAINTREP*(1.0D0/(D1*D12)-1.0D0/INTCONST)
         REPGRAD(1:3)=DUMMY*G1(1:3)
!        PRINT '(A,4I6,G15.5)','max J1,J2,REPI,REPJ,REPGRAD=',J1,J2,NREPI(J2),NREPJ(J2), &
! &                              SQRT(REPGRAD(1)**2+REPGRAD(2)**2+REPGRAD(3)**2)
         GGG(NI1+1:NI1+3)=GGG(NI1+1:NI1+3)+REPGRAD(1:3)
         GGG(NJ1+1:NJ1+3)=GGG(NJ1+1:NJ1+3)-REPGRAD(1:3)
      ENDIF
      DUMMY=0.0D0
!     IF ((.NOT.NOINT).AND.(DINT.LT.INTCONSTRAINREPCUT)) THEN
      IF ((.NOT.NOINT).AND.(DINT.LT.NREPCUT(J2))) THEN
!        D12=DSQI**6
         D12=DSQI
!        DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(12.0D0*DINT-13.0D0*INTCONSTRAINREPCUT)/INTCONST)
!        DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(12.0D0*DINT-13.0D0*NREPCUT(J2))/INTCONST)
         DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(2.0D0*DINT-3.0D0*NREPCUT(J2))/INTCONST)
         EREP=EREP+DUMMY
!        IF (DUMMY.GT.1.0D7) PRINT '(A,3I6,3G20.10)','J2,NREPI(J2),NREPJ(J2),DINT,NREPCUT(J2),DUMMY=', &
! &                                                   J2,NREPI(J2),NREPJ(J2),DINT,NREPCUT(J2),DUMMY
!        IF (((NREPI(J2).EQ.143).AND.(NREPJ(J2).EQ.191)).OR.  &
! &          ((NREPJ(J2).EQ.143).AND.(NREPI(J2).EQ.191))) THEN
!            PRINT '(A,3I6,3G20.10)','J2,NREPI(J2),NREPJ(J2),DINT,NREPCUT(J2),DUMMY=', &
! &                                   J2,NREPI(J2),NREPJ(J2),DINT,NREPCUT(J2),DUMMY
!        ENDIF
         IF (J1.EQ.2) THEN
            EEE(J1)=EEE(J1)+DUMMY
            REPEINT(J1)=REPEINT(J1)+DUMMY
            NREPINT(J1)=NREPINT(J1)+1
         ELSE IF (J1.LT.INTIMAGE+2) THEN
            EEE(J1)=EEE(J1)+DUMMY/2.0D0
            EEE(J1-1)=EEE(J1-1)+DUMMY/2.0D0
            REPEINT(J1)=REPEINT(J1)+DUMMY/2.0D0
            REPEINT(J1-1)=REPEINT(J1-1)+DUMMY/2.0D0
            NREPINT(J1)=NREPINT(J1)+1
            NREPINT(J1-1)=NREPINT(J1-1)+1
         ELSE IF (J1.EQ.INTIMAGE+2) THEN
            EEE(J1-1)=EEE(J1-1)+DUMMY
            REPEINT(J1-1)=REPEINT(J1-1)+DUMMY
            NREPINT(J1-1)=NREPINT(J1-1)+1
         ENDIF
!        DUMMY=-12.0D0*INTCONSTRAINTREP*(1.0D0/(DINT*D12)-1.0D0/INTCONST)
         DUMMY=-2.0D0*INTCONSTRAINTREP*(1.0D0/(DINT*D12)-1.0D0/INTCONST)
         REPGRAD(1:3)=DUMMY*G1INT(1:3)
!        PRINT '(A,4I6,2G15.5)','in1 J1,J2,REPI,REPJ,REPGRAD,NREPCUT=',J1,J2,NREPI(J2),NREPJ(J2), &
! &                              SQRT(REPGRAD(1)**2+REPGRAD(2)**2+REPGRAD(3)**2),NREPCUT(J2)
         GGG(NI1+1:NI1+3)=GGG(NI1+1:NI1+3)+REPGRAD(1:3)
         GGG(NJ1+1:NJ1+3)=GGG(NJ1+1:NJ1+3)-REPGRAD(1:3)
         REPGRAD(1:3)=DUMMY*G2INT(1:3)
!        PRINT '(A,4I6,2G15.5)','in1 J1,J2,REPI,REPJ,REPGRAD,NREPCUT=',J1,J2,NREPI(J2),NREPJ(J2), &
! &                              SQRT(REPGRAD(1)**2+REPGRAD(2)**2+REPGRAD(3)**2),NREPCUT(J2)
         GGG(NI2+1:NI2+3)=GGG(NI2+1:NI2+3)+REPGRAD(1:3)
         GGG(NJ2+1:NJ2+3)=GGG(NJ2+1:NJ2+3)-REPGRAD(1:3)
      ENDIF
   ENDDO
ENDDO

!
! Set gradients on frozen atoms to zero.
!
IF (FREEZE) THEN
   DO J1=2,INTIMAGE+1  
      DO J2=1,NATOMS
         IF (FROZEN(J2)) THEN
            GGG(NOPT*(J1-1)+3*(J2-1)+1)=0.0D0
            GGG(NOPT*(J1-1)+3*(J2-1)+2)=0.0D0
            GGG(NOPT*(J1-1)+3*(J2-1)+3)=0.0D0
         ENDIF
      ENDDO
   ENDDO
ENDIF
!
! Set gradients to zero for start and finish images.
!
IF (INTIMAGE.GT.0) THEN
   GGG(1:NOPT)=0.0D0
   GGG((INTIMAGE+1)*NOPT+1:(INTIMAGE+2)*NOPT)=0.0D0
ENDIF
RMS=0.0D0
RMSIMAGE(1:INTIMAGE+2)=0.0D0
DO J1=2,INTIMAGE+1
   DO J2=1,NOPT
      RMSIMAGE(J1)=RMSIMAGE(J1)+GGG(NOPT*(J1-1)+J2)**2
   ENDDO
   RMS=RMS+RMSIMAGE(J1)
   IF (LPRINT) PRINT '(A,I6,2G20.10,L5)',' congrad2> J1,EEE,RMSIMAGE,freeze=', &
  &                                                  J1,EEE(J1),RMSIMAGE(J1),IMGFREEZE(J1)
ENDDO
IF (INTIMAGE.NE.0) THEN
   RMS=SQRT(RMS/(NOPT*INTIMAGE))
ENDIF
!
! For INTIMAGE images there are INTIMAGE+2 replicas including the end points,
! and INTIMAGE+1 line segements, with associated energies stored in EEE(2:INTIMAGE+2)
!
ETOTAL=0.0D0
MAXINT=-1.0D100
MININT=1.0D100
DO J1=2,INTIMAGE+1
   ETOTAL=ETOTAL+EEE(J1)
!  IF (DEBUG) PRINT '(A,I6,A,4G15.5)',' congrad2> con/rep and con/rep int image ', &
! &      J1,' ',CONE(J1),REPE(J1),CONEINT(J1),REPEINT(J1)
   IF (CONEINT(J1)+REPEINT(J1).LT.MININT) THEN
      MININT=CONEINT(J1)+REPEINT(J1)
      NMININT=J1
   ENDIF
   IF (CONEINT(J1)+REPEINT(J1).GT.MAXINT) THEN
      MAXINT=CONEINT(J1)+REPEINT(J1)
      NMAXINT=J1
   ENDIF
ENDDO
! IF (DEBUG) PRINT '(A,G20.10,A,2I6)',' congrad2> largest  internal energy=',MAXINT,' for image ',NMAXINT
! IF (DEBUG) PRINT '(A,G20.10,A,2I6)',' congrad2> smallest internal energy=',MININT,' for image ',NMININT
IF (INTIMAGE.EQ.0) ETOTAL=EEE(1)+EEE(2)

END SUBROUTINE CONGRAD2

SUBROUTINE INTMINONLY(R1AX,R1AY,R1AZ,R2AX,R2AY,R2AZ,R1BX,R1BY,R1BZ,R2BX,R2BY,R2BZ,DINT,NOINT)
IMPLICIT NONE
DOUBLE PRECISION R1AX,R1AY,R1AZ,R2AX,R2AY,R2AZ,R1BX,R1BY,R1BZ,R2BX,R2BY,R2BZ,DINT,DUMMY
DOUBLE PRECISION DSQI, r1apr2bmr2amr1bsq, r1amr1bsq, r2amr2bsq, r1amr1bdr2amr2b, r1amr1bdr2amr2bsq
LOGICAL NOINT
!
! Is there an internal extremum?
!
! PRINT '(A,4G20.10)','r1ax,r1bx,r2ax,r2bx=',r1ax,r1bx,r2ax,r2bx
! PRINT '(A,G20.10)','(r1ax-r1bx-r2ax+r2bx)**2=',(r1ax-r1bx-r2ax+r2bx)**2
! PRINT '(A,4G20.10)','r1ay,r1by,r2ay,r2by=',r1ay,r1by,r2ay,r2by
! PRINT '(A,G20.10)','(r1ay-r1by-r2ay+r2by)**2=',(r1ay-r1by-r2ay+r2by)**2
! PRINT '(A,4G20.10)','r1az,r1bz,r2az,r2bz=',r1az,r1bz,r2az,r2bz
! PRINT '(A,G20.10)','(r1az-r1bz-r2az+r2bz)**2=',(r1az-r1bz-r2az+r2bz)**2
r1apr2bmr2amr1bsq=(r1ax-r1bx-r2ax+r2bx)**2+(r1ay-r1by-r2ay+r2by)**2+(r1az-r1bz-r2az+r2bz)**2
NOINT=.TRUE.
DINT=1.0D100
IF (r1apr2bmr2amr1bsq.EQ.0.0D0) THEN
   RETURN ! just to skip the internal solution
ELSE
   DUMMY=((r1ax-r1bx)*(r1ax-r1bx-r2ax+r2bx)+ &
 &      (r1ay-r1by)*(r1ay-r1by-r2ay+r2by)+(r1az-r1bz)*(r1az-r1bz-r2az+r2bz))/r1apr2bmr2amr1bsq
ENDIF
IF ((DUMMY.GT.0.0D0).AND.(DUMMY.LT.1.0D0)) NOINT=.FALSE.
IF (.NOT.NOINT) THEN
   r1amr1bdr2amr2b=(r1ax-r1bx)*(r2ax-r2bx)+(r1ay-r1by)*(r2ay-r2by)+(r1az-r1bz)*(r2az-r2bz)
   r1amr1bdr2amr2bsq=r1amr1bdr2amr2b**2
   r1amr1bsq=(r1ax - r1bx)**2 + (r1ay - r1by)**2 + (r1az - r1bz)**2
   r2amr2bsq=(r2ax - r2bx)**2 + (r2ay - r2by)**2 + (r2az - r2bz)**2
   DSQI=MAX((-r1amr1bdr2amr2bsq + r1amr1bsq*r2amr2bsq)/r1apr2bmr2amr1bsq,0.0D0)
   DINT=SQRT(DSQI)
ENDIF

RETURN

END SUBROUTINE INTMINONLY
