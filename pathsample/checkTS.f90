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
! Using input arguments to determine whether we should discount this TS or not.

      SUBROUTINE CHECKTS(LETS,LEPLUS,LEMINUS,LKPLUS,LKMINUS,LNCONNPLUS,LNCONNMINUS,LPLUS,LMINUS,DEGENT,CUT_UNDERFLOW,DEADTS)
      USE PORFUNCS
      USE COMMON, ONLY : NCONNMIN, TSTHRESH, MAXBARRIER
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: LNCONNPLUS, LNCONNMINUS, LPLUS, LMINUS
      DOUBLE PRECISION, INTENT(IN) :: LETS, LEPLUS, LEMINUS, LKPLUS, LKMINUS, CUT_UNDERFLOW
      LOGICAL, INTENT(IN) :: DEGENT
      LOGICAL, INTENT(OUT) :: DEADTS

      IF ((LNCONNPLUS.LE.NCONNMIN).OR.(LNCONNMINUS.LE.NCONNMIN).OR.(LETS.GT.TSTHRESH).OR. &
          (LKPLUS.LT.CUT_UNDERFLOW).OR.(LKMINUS.LT.CUT_UNDERFLOW)) THEN
         DEADTS=.TRUE.
      ELSEIF ((LETS-LEPLUS.GT.MAXBARRIER).AND.(LETS-LEMINUS.GT.MAXBARRIER)) THEN
         DEADTS=.TRUE.
      ELSEIF ((LPLUS.EQ.LMINUS).AND.(.NOT.DEGENT)) THEN
         DEADTS=.TRUE.
      ELSE 
         DEADTS=.FALSE.
      ENDIF

      RETURN
      END SUBROUTINE CHECKTS
