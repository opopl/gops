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
C
C the function DF1DIM has to accompany LINMIN
C
      SUBROUTINE ZWISCHEN(XX,DF1DIM,POTEL1)
      USE commons
      use f1com
      IMPLICIT DOUBLE PRECISION (A-H,P-Z)
      DOUBLE PRECISION XT(3*NATOMS),DF(3*NATOMS)

      DO 10 J=1,NCOM
         XT(J)=PCOM(J)+XX*XICOM(J)
10    CONTINUE
      CALL POTENTIAL(XT,DF,POTEL,.TRUE.,.FALSE.)
      POTEL1=POTEL
      DF1DIM=0.0D0
      DO 20 J=1,NCOM
         DF1DIM=DF1DIM+DF(J)*XICOM(J)
20    CONTINUE
C     PRINT*,'energy in zwischen=',POTEL
      RETURN
      END
