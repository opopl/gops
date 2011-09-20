!op226>=================================== 
!op226> GPL License Info {{{ 
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
!op226>}}} 
!op226>=================================== 
!> \brief This subroutine ...
!> \name MCRUNS
!> \param SCREENC  
!op226>=================================== 
      SUBROUTINE MCRUNS(SCREENC)
      USE COMMONS
      IMPLICIT NONE
      LOGICAL LOPEN
      DOUBLE PRECISION SCREENC(3*NATOMS)

      INTEGER J1

      INQUIRE(UNIT=1,OPENED=LOPEN)

      IF (LOPEN) THEN
         WRITE(*,'(A,I2,A)') 'mcruns> ERROR *** Unit ', 1, ' is not free '
         STOP
      ENDIF

      DO J1=1,NRUNS
         CALL MC(MCSTEPS(J1),TFAC(J1),SCREENC)
      ENDDO

      CALL FINALQ
      CALL FINALIO

      RETURN
      END
