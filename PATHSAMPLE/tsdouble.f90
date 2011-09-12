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

SUBROUTINE TSDOUBLE
   USE COMMON
   IMPLICIT NONE
   DOUBLE PRECISION, ALLOCATABLE :: VDP(:)
   INTEGER, ALLOCATABLE :: VINT(:)

   ALLOCATE(VDP(MAXTS),VINT(MAXTS))

   VDP(1:MAXTS)=ETS(1:MAXTS)
   DEALLOCATE(ETS)
   ALLOCATE(ETS(2*MAXTS))
   ETS(1:MAXTS)=VDP(1:MAXTS)

   VDP(1:MAXTS)=FVIBTS(1:MAXTS)
   DEALLOCATE(FVIBTS)
   ALLOCATE(FVIBTS(2*MAXTS))
   FVIBTS(1:MAXTS)=VDP(1:MAXTS)

   VDP(1:MAXTS)=KPLUS(1:MAXTS)
   DEALLOCATE(KPLUS)
   ALLOCATE(KPLUS(2*MAXTS))
   KPLUS(1:MAXTS)=VDP(1:MAXTS)

   VDP(1:MAXTS)=KMINUS(1:MAXTS)
   DEALLOCATE(KMINUS)
   ALLOCATE(KMINUS(2*MAXTS))
   KMINUS(1:MAXTS)=VDP(1:MAXTS)

   VDP(1:MAXTS)=IXTS(1:MAXTS)
   DEALLOCATE(IXTS)
   ALLOCATE(IXTS(2*MAXTS))
   IXTS(1:MAXTS)=VDP(1:MAXTS)

   VDP(1:MAXTS)=IYTS(1:MAXTS)
   DEALLOCATE(IYTS)
   ALLOCATE(IYTS(2*MAXTS))
   IYTS(1:MAXTS)=VDP(1:MAXTS)

   VDP(1:MAXTS)=IZTS(1:MAXTS)
   DEALLOCATE(IZTS)
   ALLOCATE(IZTS(2*MAXTS))
   IZTS(1:MAXTS)=VDP(1:MAXTS)

   VDP(1:MAXTS)=NEGEIG(1:MAXTS)
   DEALLOCATE(NEGEIG)
   ALLOCATE(NEGEIG(2*MAXTS))
   NEGEIG(1:MAXTS)=VDP(1:MAXTS)

   VINT(1:MAXTS)=HORDERTS(1:MAXTS)
   DEALLOCATE(HORDERTS)
   ALLOCATE(HORDERTS(2*MAXTS))
   HORDERTS(1:MAXTS)=VINT(1:MAXTS)

   VINT(1:MAXTS)=PLUS(1:MAXTS)
   DEALLOCATE(PLUS)
   ALLOCATE(PLUS(2*MAXTS))
   PLUS(1:MAXTS)=VINT(1:MAXTS)

   VINT(1:MAXTS)=MINUS(1:MAXTS)
   DEALLOCATE(MINUS)
   ALLOCATE(MINUS(2*MAXTS))
   MINUS(1:MAXTS)=VINT(1:MAXTS)

   VINT(1:MAXTS)=POINTERM(1:MAXTS)
   DEALLOCATE(POINTERM)
   ALLOCATE(POINTERM(2*MAXTS))
   POINTERM(1:MAXTS)=VINT(1:MAXTS)

   VINT(1:MAXTS)=POINTERP(1:MAXTS)
   DEALLOCATE(POINTERP)
   ALLOCATE(POINTERP(2*MAXTS))
   POINTERP(1:MAXTS)=VINT(1:MAXTS)

   IF (DIJKSTRAT .OR. KSHORTESTPATHST) THEN
      VINT(1:MAXTS)=TSATTEMPT(1:MAXTS)
      DEALLOCATE(TSATTEMPT)
      ALLOCATE(TSATTEMPT(2*MAXTS))
      TSATTEMPT(1:MAXTS)=VINT(1:MAXTS)
   ENDIF

   MAXTS=2*MAXTS
   PRINT '(A,I8)','tsdouble> Maximum number of ts increased to ',MAXTS
   DEALLOCATE(VDP,VINT)

END SUBROUTINE TSDOUBLE
