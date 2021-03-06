!   OPTIM: A program for optimizing geometries and calculating reaction pathways
!   Copyright (C) 1999-2006 David J. Wales
!   This file is part of OPTIM.
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
SUBROUTINE MYCPU_TIME(POO,CHECKLIMIT)
use porfuncs
USE KEY
IMPLICIT NONE
DOUBLE PRECISION MYTIME
DOUBLE PRECISION, INTENT(OUT) :: POO
LOGICAL, INTENT(IN) :: CHECKLIMIT

CALL CPU_TIME(MYTIME)
POO=MYTIME ! without this extra assignment NAG f95 returns a random number!
           ! saving TSTART is necessary for PG compiler, where initial time is
           ! not zero!
IF (.NOT.CHECKLIMIT) RETURN ! we don;t want to stop half way through a path, for example
IF (MYTIME-TSTART.GT.TIMELIMIT) THEN
   PRINT '(A,F12.1,A,F12.1)','CPU time limit exceeded: time=',MYTIME-TSTART,' limit=',TIMELIMIT
   IF (DUMPSP) CALL DODUMPSP
   STOP
ENDIF

RETURN
END SUBROUTINE MYCPU_TIME
