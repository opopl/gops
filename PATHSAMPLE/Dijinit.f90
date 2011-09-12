!   PATHSAMPLE: A driver for OPTIM to create stationary point databases using discrete path sampling and perform kinetic analysis
!   Copyright (C) 1999-2009 David J. Wales
!   This file is part of PATHSAMPLE.
!
!   PATHSAMPLE is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   PATHSAMPLE is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  Dijkstra connection algorithm for pathsample.
!
SUBROUTINE DIJINIT(NWORST)
USE PORFUNCS
USE COMMON
IMPLICIT NONE

INTEGER J1, J2, J4, PARENT(NMIN), JMINW, NPERM, J5, LJ1, LJ2, NWORST, NSTEPS, NMINSTART, NMINEND
INTEGER, ALLOCATABLE :: LOCATIONSTART(:), LOCATIONEND(:)
LOGICAL PERMANENT(NMIN), ISA(NMIN), ISB(NMIN), ISSTART(NMIN), NOTDONE
DOUBLE PRECISION MINWEIGHT, TMPWEIGHT, WEIGHT(NMIN), DUMMY, TNEW, ELAPSED, PFTOTALSTART, HUGESAVE, THRESH
DOUBLE PRECISION MAXWEIGHT, SCALEFAC

CALL CPU_TIME(ELAPSED)
DEALLOCATE(DMIN1,DMIN2)
ALLOCATE(DMIN1(10000),DMIN2(10000)) ! surely this will be enough room!
ISA(1:NMIN)=.FALSE.
ISB(1:NMIN)=.FALSE.
DO J1=1,NMINA
   ISA(LOCATIONA(J1))=.TRUE.
ENDDO
DO J1=1,NMINB
   ISB(LOCATIONB(J1))=.TRUE.
ENDDO

!!!!!!!!!!!!!!!!!!!   Dijkstra calculation    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Dijkstra connect process similar to that in OPTIM with a weight for missing
!  connections based on a distance metric.
!
IF (DIRECTION.EQ.'AB') THEN
   NMINSTART=NMINB
   NMINEND=NMINA
   ALLOCATE(LOCATIONSTART(NMINB),LOCATIONEND(NMINA))
   LOCATIONSTART(1:NMINB)=LOCATIONB(1:NMINB)
   LOCATIONEND(1:NMINA)=LOCATIONA(1:NMINA)
   ISSTART(1:NMIN)=ISB(1:NMIN)
   PFTOTALSTART=PFTOTALB
ELSEIF (DIRECTION.EQ.'BA') THEN
   NMINSTART=NMINA
   NMINEND=NMINB
   ALLOCATE(LOCATIONSTART(NMINA),LOCATIONEND(NMINB))
   LOCATIONSTART(1:NMINA)=LOCATIONA(1:NMINA)
   LOCATIONEND(1:NMINB)=LOCATIONB(1:NMINB)
   ISSTART(1:NMIN)=ISA(1:NMIN)
   PFTOTALSTART=PFTOTALA
ENDIF
!
!  Find largest weight for each B(A) minimum to all A(B) minima.
!
!  Added maximum weight condition via a distance scaling. 
!  Otherwise loss of precision can cause connections to be missed completely. DJW 29/7/08
!
loopstart: DO J1=1,NMINSTART ! cycle over all minima in the starting state
   MAXWEIGHT=1.0D10
   SCALEFAC=1.0D0
222   LJ1=LOCATIONSTART(J1)
   WEIGHT(1:NMIN)=HUGE(1.0D0)
   HUGESAVE=WEIGHT(1)
   WEIGHT(LJ1)=0.0D0
   PERMANENT(1:NMIN)=.FALSE.
   PERMANENT(LJ1)=.TRUE.
   NPERM=1
   PARENT(1:NMIN)=0 ! parent is initially undefined
   J4=LJ1

   dijkstraloop: DO

      DO J2=1,NMIN
         IF (J2.EQ.J4) CYCLE
         IF (PERMANENT(J2)) CYCLE
         TMPWEIGHT=PAIRDIST(MAX(J2,J4)*(MAX(J2,J4)-1)/2+MIN(J4,J2))*SCALEFAC
         IF (TMPWEIGHT.LT.HUGE(1.0D0)/10.0D0) THEN ! don;t raise a huge number to any power!
            IF (INDEXCOSTFUNCTION) THEN 
               IF (TMPWEIGHT.EQ.0.0D0) THEN ! minima are connected!
               ELSE
                  TMPWEIGHT=ABS(J4-J2)
                  IF (DIRECTION.EQ.'BA') THEN
                     IF (J4.LE.NMINA) TMPWEIGHT=NMIN+1-J2 ! not sure that this really makes sense for A and B ! DJW
                     IF (J2.LE.NMINA) TMPWEIGHT=NMIN+1-J4
                  ELSE
                     IF ((J4.LE.NMINA+NMINB).AND.(J4.GT.NMINA)) TMPWEIGHT=NMIN+1-J2
                     IF ((J2.LE.NMINA+NMINB).AND.(J2.GT.NMINA)) TMPWEIGHT=NMIN+1-J4
                  ENDIF
                ENDIF
            ELSEIF (EXPCOSTFUNCTION) THEN ! saves memory and CPU when endpoint separation is very large SAT
               IF (TMPWEIGHT.EQ.0.0D0) THEN
                  ! do nothing - don;t set the weight to one !! DJW 22/7/08
               ELSEIF (TMPWEIGHT.GT.300.0D0) THEN
                  TMPWEIGHT=EXP(300.0D0)
               ELSE
                  TMPWEIGHT=EXP(TMPWEIGHT)
               ENDIF
            ELSE ! compare squares to favour more small jumps over big ones DJW
               IF (TMPWEIGHT.EQ.0.0D0) THEN
               ELSEIF (COSTFUNCTIONPOWER.EQ.0) THEN
                  TMPWEIGHT=TMPWEIGHT+1.0D0
               ELSEIF (COSTFUNCTIONPOWER.EQ.-1) THEN
                  TMPWEIGHT=1.0D0/TMPWEIGHT
               ELSE
                  TMPWEIGHT=TMPWEIGHT**COSTFUNCTIONPOWER 
               ENDIF
            ENDIF
         ENDIF
         
         IF (TMPWEIGHT+WEIGHT(J4).LT.WEIGHT(J2)) THEN ! relax J2
            WEIGHT(J2)=TMPWEIGHT+WEIGHT(J4)
            PARENT(J2)=J4
         ENDIF
      ENDDO

      MINWEIGHT=HUGE(1.0D0)
      NOTDONE=.TRUE.
      DO J2=1,NMIN
         IF (.NOT.PERMANENT(J2)) THEN
            IF (WEIGHT(J2).LT.MINWEIGHT) THEN
               MINWEIGHT=WEIGHT(J2)
               JMINW=J2
               NOTDONE=.FALSE.
            ENDIF
         ENDIF
      ENDDO
      IF (NOTDONE) PRINT '(A,I8,A,I8)','dijinit> WARNING - JMINW not set - value=',JMINW,' J4=',J4

      J4=JMINW
      PERMANENT(J4)=.TRUE.
      NPERM=NPERM+1
      IF (WEIGHT(J4).GT.MAXWEIGHT) THEN
         SCALEFAC=SCALEFAC/2.0D0
         PRINT '(A,G20.10)','dijinit> Maximum weight is too large - scaling distances by ',SCALEFAC
         GOTO 222
      ENDIF

      IF (NPERM.EQ.NMIN) EXIT dijkstraloop

   ENDDO dijkstraloop

ENDDO loopstart
! 
!  Summarise the best path for any A(B) and any B(A)
!
LJ2=LOCATIONEND(1)
LJ1=LOCATIONSTART(1)
J5=LJ2
NWORST=0
NSTEPS=0
PRINT '(A)','Dijinit> Summary of best path based on missing connection metric - note distance scaling is removed'
PRINT '(A)','    min1          energy        min2          energy             metric          edge weight            weight'
DO 
   IF (PARENT(J5).EQ.0) THEN
      PRINT '(A,I6,A)','Dijinit> ERROR - parent for J5=',J5,' is zero'
      PRINT '(A)',     'Dijinit> Suggests all possible pairs have been tried!'
      STOP
   ENDIF
   DUMMY=PAIRDIST(MAX(J5,PARENT(J5))*(MAX(J5,PARENT(J5))-1)/2+MIN(J5,PARENT(J5)))*SCALEFAC
   IF (DUMMY.LT.HUGE(1.0D0)/10.0D0) THEN ! don;t raise a huge number to any power!
      IF (INDEXCOSTFUNCTION) THEN
         IF (DUMMY.EQ.0.0D0) THEN ! minima are connected!
            TMPWEIGHT=0.0D0
         ELSE
            TMPWEIGHT=ABS(J5-PARENT(J5))
            IF (DIRECTION.EQ.'AB') THEN
               IF (J5.LE.NMINA) TMPWEIGHT=NMIN+1-PARENT(J5)
               IF (PARENT(J5).LE.NMINA) TMPWEIGHT=NMIN+1-J5
            ELSE
               IF ((PARENT(J5).LE.NMINA+NMINB).AND.(PARENT(J5).GT.NMINA)) TMPWEIGHT=NMIN+1-J5
               IF ((J5.LE.NMINA+NMINB).AND.(J5.GT.NMINA)) TMPWEIGHT=NMIN+1-PARENT(J5)
            ENDIF 
          ENDIF
      ELSEIF (EXPCOSTFUNCTION) THEN ! saves memory and CPU when endpoint separation is very large SAT
         IF (DUMMY.EQ.0.0D0) THEN
            TMPWEIGHT=0.0D0
         ELSEIF (DUMMY.GT.300.0D0) THEN
             TMPWEIGHT=EXP(300.0D0)
         ELSE
            TMPWEIGHT=EXP(DUMMY)
         ENDIF
      ELSE ! compare higher powers to favour more small jumps over big ones DJW
         IF (DUMMY.EQ.0.0D0) THEN
            TMPWEIGHT=0.0D0
         ELSEIF (INTERPCOSTFUNCTION) THEN
            TMPWEIGHT=DUMMY**COSTFUNCTIONPOWER
         ELSEIF (COSTFUNCTIONPOWER.EQ.0) THEN
            TMPWEIGHT=1.0D0
         ELSEIF (COSTFUNCTIONPOWER.EQ.-1) THEN
            TMPWEIGHT=1.0D0/TMPWEIGHT
         ELSE
            TMPWEIGHT=DUMMY**COSTFUNCTIONPOWER
         ENDIF
      ENDIF
   ELSE
      TMPWEIGHT=DUMMY
   ENDIF
   NSTEPS=NSTEPS+1
   
   PRINT '(2(I8,G20.10),3G20.10)',J5,EMIN(J5),parent(J5),EMIN(PARENT(J5)),DUMMY/SCALEFAC,TMPWEIGHT,WEIGHT(J5)
   THRESH=0.0D0
   IF (BHINTERPT) THRESH=BHDISTTHRESH ! for bhinterp runs raise the threshold to BHDISTTHRESH
   IF (BISECTT) THRESH=BISECTMINDIST ! for bisect runs raise the threshold to BISECTMINDIST
   IF ((DUMMY/SCALEFAC.GT.THRESH).AND.(TMPWEIGHT.LT.HUGE(1.0D0)/10.0D0)) THEN
      NWORST=NWORST+1
      IF (NWORST.GT.10000) THEN
         PRINT '(A,I8)','ERROR in Dijinit, too many gaps, NWORST=',NWORST
         STOP
      ENDIF
      DMIN1(NWORST)=J5
      DMIN2(NWORST)=PARENT(J5)
!
!  On the first attempt there may only be two minima, and with more than one cpu available
!  we may need to choose the same best path several times before we can get
!  more than one distinct candidate.
!
      IF (NMIN.GT.2) PAIRDIST(MAX(J5,PARENT(J5))*(MAX(J5,PARENT(J5))-1)/2+MIN(J5,PARENT(J5)))=HUGE(1.0D0)
   ENDIF
   J5=PARENT(J5)
   IF (J5.EQ.LJ1) EXIT
   IF (J5.EQ.0) EXIT
ENDDO
PRINT '(2(A,I8))','Dijinit> Number of steps=',NSTEPS,' number of missing connections=',NWORST
IF (NWORST.EQ.0) THEN
   PRINT '(A)','Dijinit> Connected path found'
!  IF (DIJCONT) THEN
!     DIJINITT=.FALSE.
!  ELSE
      STOP
!  ENDIF
ENDIF
CALL CPU_TIME(TNEW)
TDIJKSTRA=TDIJKSTRA+TNEW-ELAPSED

RETURN
END
