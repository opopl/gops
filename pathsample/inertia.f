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

! ################################################################################

! Wrapper for inertia calculation to avoid explicit coordinate transformations
! where possible

      subroutine inertiaWrapper(Q, nCoordSites, angleAxis, ITX, ITY, ITZ)

      use rigidBodymod

      implicit none

! Subroutine arguments

      integer, intent(IN) :: nCoordSites
      real (kind=kind(0.0d0)), intent(INOUT) :: Q(3*nCoordSites)
      logical, intent(IN) :: angleAxis
      real (kind=kind(0.0d0)), intent(OUT) :: ITX, ITY, ITZ

! Local variables

      integer :: numCartPoints
      real (kind=kind(0.0d0)), allocatable :: cartCoords(:)

      if (angleAxis) then
         numCartPoints = (nCoordSites/2)*rbPotential%nPhysicalSites
         allocate (cartCoords(numCartPoints*3))
         call systemToCartesians(nCoordSites/2, Q, cartCoords)
         call inertia(cartCoords,numCartPoints,ITX,ITY,ITZ)
         deallocate (cartCoords)
      else
         call inertia(Q,nCoordSites,ITX,ITY,ITZ)
      endif
!     PRINT *,'ccords in inertiawrapper:'
!     PRINT '(3F20.10)',Q(1:3*nCoordSites)
!     PRINT '(A,3F20.10)','moments of inertia: ',ITX,ITY,ITZ

      end

C*************************************************************************
C
C  Build inertia tensor and get principal values.
C
      SUBROUTINE INERTIA(Q,nCartPoints,ITX,ITY,ITZ)
      USE COMMON
      IMPLICIT NONE
      integer, intent(IN) :: nCartPoints
      INTEGER J1, J2, J3
      DOUBLE PRECISION IT(3,3), Q(3*nCartPoints), CMX, CMY, CMZ, VEC(3,3), ITX, ITY, ITZ, MASST

      if (size(MASS).ne.nCartPoints) then
         print *, 'inertia> Size of MASS not equal to number of points'
         stop
      endif

      IF (RBAAT) NATOMS = NATOMS/2

      CMX=0.0D0
      CMY=0.0D0
      CMZ=0.0D0
      MASST=0.0D0
      DO J1=1,nCartPoints
         CMX=CMX+Q(3*(J1-1)+1)*MASS(J1)
         CMY=CMY+Q(3*(J1-1)+2)*MASS(J1)
         CMZ=CMZ+Q(3*(J1-1)+3)*MASS(J1)
         MASST=MASST+MASS(J1)
C        WRITE(*,'(A,I5,4F20.10)') 'I,ATMASS(I),=',J1,MASS(J1),Q(3*(J1-1)+1),Q(3*(J1-1)+2),Q(3*(J1-1)+3)
      ENDDO
      CMX=CMX/MASST
      CMY=CMY/MASST
      CMZ=CMZ/MASST
C     PRINT*,'MASST,CMX,CMY,CMZ=',MASST,CMX,CMY,CMZ
      DO J1=1,nCartPoints
         Q(3*(J1-1)+1)=Q(3*(J1-1)+1)-CMX
         Q(3*(J1-1)+2)=Q(3*(J1-1)+2)-CMY
         Q(3*(J1-1)+3)=Q(3*(J1-1)+3)-CMZ
      ENDDO

      DO J1=1,3
         DO J2=1,3
            IT(J1,J2)=0.0D0
            DO J3=1,nCartPoints
               IT(J1,J2)=IT(J1,J2)-Q(3*(J3-1)+J1)*Q(3*(J3-1)+J2)*MASS(J3)
            ENDDO
            IF (J1.EQ.J2) THEN
               DO J3=1,nCartPoints
                  IT(J1,J2)=IT(J1,J2)+(Q(3*(J3-1)+1)**2+Q(3*(J3-1)+2)**2+Q(3*(J3-1)+3)**2)*MASS(J3)
               ENDDO
            ENDIF
         ENDDO
      ENDDO

      CALL EIG(IT,VEC,3,3,0)
      ITX=IT(1,1)
      ITY=IT(2,2)
      ITZ=IT(3,3)
!
!  Must not move frozen atoms!
!
      DO J1=1,nCartPoints
         Q(3*(J1-1)+1)=Q(3*(J1-1)+1)+CMX
         Q(3*(J1-1)+2)=Q(3*(J1-1)+2)+CMY
         Q(3*(J1-1)+3)=Q(3*(J1-1)+3)+CMZ
      ENDDO

      IF (RBAAT) NATOMS = 2*NATOMS

      RETURN
      END
