!   GMIN: A program for finding global minima
!   Copyright (C) 1999-2010 David J. Wales
!   This file is part of GMIN.
!
!   GMIN is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   GMIN is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
! 
! The GROUPROTATION subroutine allows for almost any rotation of a defined set of atoms.
! The rotation axis is defined by two atoms (BATOMS1 and BATOM2), the group of atoms to 
! rotate is defined by the logical array ATOMINGROUP (if element is .TRUE., that atom is
! in the group), ANGLE is the rotation angle in radians and STEPCOORDS contains the current
! atomic coordinates.
!

! NEW ROUTINE TO DRIVE GROUPROTATION MOVES
! JP is the parallel run ID so that only the appropriate coordinates are
! altered during parallel runs.
      SUBROUTINE GROUPROTSTEP(JP)
      USE commons
      IMPLICIT NONE
      DOUBLE PRECISION :: DPRAND, PI, TWOPI, GROUPROTANGLE, GROUPROTANGLEDEG
      INTEGER :: I1,JP
! Some helpful parameters
      PI=ATAN(1.0D0)*4
      TWOPI=2.0D0*PI
! For each group....      
      DO I1=1,NGROUPS
         IF(ATOMGROUPPSELECT(I1).GE.DPRAND()) THEN
! Group selected to be rotated - calculate rotation angle
            GROUPROTANGLE=(DPRAND()-0.5)*twopi*ATOMGROUPSCALING(I1)
            GROUPROTANGLEDEG=GROUPROTANGLE*(180/pi)
! Print some into to GMIN_out for the user
            WRITE(LFH,*) 'GROUPROTATION> Rotating group ',TRIM(ADJUSTL(ATOMGROUPNAMES(I1))),' by ',GROUPROTANGLEDEG
! Call the rotation subroutine
            CALL GROUPROTATION(ATOMGROUPAXIS(I1,1),ATOMGROUPAXIS(I1,2),GROUPROTANGLE,ATOMGROUPS(I1,:),COORDS(:,JP))
         ENDIF
      ENDDO 
      END SUBROUTINE GROUPROTSTEP

      SUBROUTINE GROUPROTATION(BATOM1,BATOM2,ANGLE,ATOMINGROUP,STEPCOORDS)
      USE commons
      IMPLICIT NONE
      INTEGER :: BATOM1, BATOM2, I1
      DOUBLE PRECISION :: BVECTOR(3), LENGTH, ANGLE, DUMMYMAT(3,3)=0.0D0, ROTMAT(3,3)
      DOUBLE PRECISION :: GROUPATOM(3), GROUPATOMROT(3), STEPCOORDS(3*NATOMS)
      LOGICAL :: ATOMINGROUP(NATOMS)
! STEP 1
! Produce notmalised bond vector corresponding to the rotation axis
! BATOM1 and BATOM2 are the atoms defining this vector
      BVECTOR(1)=STEPCOORDS(3*BATOM1-2)-STEPCOORDS(3*BATOM2-2)
      BVECTOR(2)=STEPCOORDS(3*BATOM1-1)-STEPCOORDS(3*BATOM2-1)
      BVECTOR(3)=STEPCOORDS(3*BATOM1  )-STEPCOORDS(3*BATOM2  )
! Find length   
      LENGTH=DSQRT(BVECTOR(1)**2 + BVECTOR(2)**2 + BVECTOR(3)**2)
! Normalise      
      BVECTOR(1)=BVECTOR(1)/LENGTH
      BVECTOR(2)=BVECTOR(2)/LENGTH
      BVECTOR(3)=BVECTOR(3)/LENGTH
! STEP 2
! Scale this vector so its length is the rotation to be done (in radians)
      BVECTOR(1)=BVECTOR(1)*ANGLE
      BVECTOR(2)=BVECTOR(2)*ANGLE
      BVECTOR(3)=BVECTOR(3)*ANGLE
! STEP 3
! Get the rotation matrix for this vector axis from RMDRVT
! Interface:
! SUBROUTINE RMDRVT(P, RM, DRM1, DRM2, DRM3, GTEST)
! P is an un-normalised vector you wish to rotate around. Its length equals the desired rotation in radians
! RM will return the 3x3 rotation matrix
! DRM1-3 are derivative matricies, not needed here
! GTEST is also not needed so set to .FALSE.
      CALL RMDRVT(BVECTOR,ROTMAT,DUMMYMAT,DUMMYMAT,DUMMYMAT,.FALSE.)
! STEP 4
! Rotate group, one atom at a time. First, translate atom so the pivot (end of bond closest to atom) is at the origin  
      DO I1=1,NATOMS
         IF (ATOMINGROUP(I1)) THEN
            GROUPATOM(1)=STEPCOORDS(3*I1-2)-STEPCOORDS(3*BATOM2-2)
            GROUPATOM(2)=STEPCOORDS(3*I1-1)-STEPCOORDS(3*BATOM2-1)
            GROUPATOM(3)=STEPCOORDS(3*I1  )-STEPCOORDS(3*BATOM2  )
! Apply the rotation matrix
            GROUPATOMROT=MATMUL(ROTMAT,GROUPATOM)
! Translate back to the origin and copy to COORDS
            STEPCOORDS(3*I1-2)=GROUPATOMROT(1)+STEPCOORDS(3*BATOM2-2)
            STEPCOORDS(3*I1-1)=GROUPATOMROT(2)+STEPCOORDS(3*BATOM2-1)
            STEPCOORDS(3*I1  )=GROUPATOMROT(3)+STEPCOORDS(3*BATOM2  )
         ENDIF
      ENDDO

      END SUBROUTINE GROUPROTATION
