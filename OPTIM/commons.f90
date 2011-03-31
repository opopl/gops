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
MODULE COMMONS
      IMPLICIT NONE

      INTEGER :: NATOMS, ISTATUS, NUNIQUE, NOPT, IPRNT, INR, IVEC, IVEC2, ISTCRT, NMOL, NINTS, NRBSITES
      DOUBLE PRECISION :: ZSTAR, RHO, PARAM1=0.0D0, PARAM2, PARAM3, PARAM4, PARAM5, PARAM6, PARAM7, EVDISTTHRESH, NEWRES_TEMP
      LOGICAL REDOPATH, REDOPATHXYZ, RATIOS, REDOPATHNEB

      INTEGER, ALLOCATABLE, DIMENSION(:) :: NR, IATNUM      !   MXATMS
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: ATMASS !   MXATMS
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: STPMAX ! 3*MXATMS
      CHARACTER(LEN=5), ALLOCATABLE, DIMENSION(:) :: ZSYM   !   MXATMS
      DOUBLE PRECISION, ALLOCATABLE :: TAGFAC(:)
      DOUBLE PRECISION, ALLOCATABLE :: RBSITE(:,:), RBSTLA(:,:), STCHRG(:), DPMU(:)
      INTEGER, ALLOCATABLE :: TAGNUM(:)
      INTEGER NTAG
      LOGICAL TAGT, CHANGE_TEMP

END MODULE COMMONS
