!     Copyright (C) 1999-2008 David J. Wales
!  This file is part of OPTIM.
!
!  OPTIM is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  OPTIM is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  COORDSA becomes the optimal alignment of the optimal permutation(-inversion)
!  isomer, but without the permutations. DISTANCE is the residual square distance
!  for the best alignment with respect to permutation(-inversion)s as well as
!  orientation and centre of mass.
!
!  MYORIENT is called first for both COORDSA and COORDSB to put them into
!  a standard orientation in DUMMYA and DUMMYB (which both have the centre of
!  coordinates at the origin). 
!  The objective is to identify permutation-inversion isomers without fail. 
!  However, we have to cycle over all equivalent atoms in two particular orbits for DUMMYA
!  to achieve this.
!  We iterate permutations and newmindist minimisations up to a maximum number or
!  until no more permutations are required for each instance of DUMMYA aligned 
!  according to NCHOOSE1 and NCHOOSE2 by MYORIENT. The cumulative rotation
!  matrix that takes the initial DUMMYA to the one that aligns best with DUMMYB
!  is saved in RMATCUMUL.
!  Then, if we've not going BULK, AMBER, or CHARMM, we try again for the inverted
!  version of COORDSA. The transformation corresponding to the minimum distance
!  is saved whenever it is improved - the best alignment including permutations
!  is saved in XBEST, and the last step is to rotate this back to coincide best
!  with COORDSB (rather than DUMMYB) using ROTINVBBEST. This gives suitable
!  fixed end points for DNEB.
!  Finally, we transform COORDSA to be in optimal alignment, but without the
!  permutations in XBEST. The overall transformation is
!  COORDSA -> +/- ROTINVB RMATCUMUL ROTA (COORDSA - CMA) 
!
!  The correspondence between COORDSA and DUMMYA after DUMMYA has been aligned by
!  newmindist is
!  +/- RMATCUMUL ROTA (COORDSA - CMA) = permutation(DUMMYA)
!  where +/- is given by the value of INVERT.
!  The centres of coordinates for COORDSA and COORDSB can be anywhere. On return, the
!  centre of coordinates of COORDSA will be the same as for COORDSB.
!
SUBROUTINE MINPERMDIST(COORDSB,COORDSA,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGID,RMATBEST)
USE COMMONS,ONLY : NPERMGROUP, NPERMSIZE, PERMGROUP, NSETS, SETS, CHRMMT, MYUNIT, STOCKT, NFREEZE, &
  & AMBERT, CSMGPINDEX, CSMT, ZSYM, PERIODIC, PERMDIST, PULLT, EFIELDT
USE PORFUNCS
IMPLICIT NONE

INTEGER, PARAMETER :: MAXIMUMTRIES=20
INTEGER NATOMS, NPERM, PATOMS, NTRIES, ISTAT
INTEGER J3, INVERT, NORBIT1, NORBIT2, NCHOOSE2, NDUMMY, LPERM(NATOMS), J1, J2, NCHOOSE1
DOUBLE PRECISION DIST2, COORDSA(3*NATOMS), COORDSB(3*NATOMS), DISTANCE, DUMMYA(3*NATOMS), DUMMYB(3*NATOMS), DUMMY(3*NATOMS)
DOUBLE PRECISION BOXLX,BOXLY,BOXLZ,WORSTRAD,RMAT(3,3),ENERGY, VNEW(3*NATOMS), DX, DY, DZ, RMS, DBEST, XBEST(3*NATOMS)
DOUBLE PRECISION CMXA, CMXB, CMXC, PDISTANCE
DOUBLE PRECISION ROTA(3,3), ROTINVA(3,3), ROTB(3,3), ROTINVB(3,3), ROTINVBBEST(3,3), ROTABEST(3,3), RMATBEST(3,3), TMAT(3,3)
DOUBLE PRECISION CMAX, CMAY, CMAZ, CMBX, CMBY, CMBZ, RMATCUMUL(3,3)
DOUBLE PRECISION REFXZ(3,3)
LOGICAL DEBUG, TWOD, RIGID, BULKT
DOUBLE PRECISION PDUMMYA(3*NATOMS), PDUMMYB(3*NATOMS), LDISTANCE, DUMMYC(3*NATOMS), XDUMMY, GEOMDIFFTOL
SAVE NORBIT1, NORBIT2
INTEGER NEWPERM(NATOMS), ALLPERM(NATOMS), SAVEPERM(NATOMS)

! DUMMYB(1:3*NATOMS)=COORDSB(1:3*NATOMS)
! DUMMYA(1:3*NATOMS)=COORDSA(1:3*NATOMS)

GEOMDIFFTOL=0.5D0
DISTANCE = 0.0D0
IF (.NOT.PERMDIST) THEN
   CALL NEWMINDIST(COORDSB,COORDSA,NATOMS,DISTANCE,PERIODIC,TWOD,ZSYM(1),.FALSE.,RIGID,DEBUG,RMAT)
   RETURN
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALL OCHARMM(DUMMYA,VNEW,ENERGY,.FALSE.,.FALSE.)
! CALL POTENTIAL(DUMMYA,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
! PRINT '(2(A,F25.15))',' Initial energy=',ENERGY,' RMS=',RMS
! PRINT '(2(A,F25.15))',' for coordinates:'
! PRINT '(3F25.15)',DUMMYA(1:3*NATOMS)
! PRINT '(A,F25.15,A)',' Initial energy=',ENERGY,' kcal/mol'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  OPEN(UNIT=10,FILE='minpermdist.xyz',STATUS='UNKNOWN')
!  WRITE(10,'(I6)') NATOMS/2
!  WRITE(10,'(A)') 'A initial'
!  DO J3=1,NATOMS/2
!      WRITE(10,'(A2,2X,3F20.10)') 'LA',COORDSA(3*(J3-1)+1),COORDSA(3*(J3-1)+2),COORDSA(3*(J3-1)+3)
!  ENDDO
!  WRITE(10,'(I6)') NATOMS/2
!  WRITE(10,'(A)') 'B initial'
!  DO J3=1,NATOMS/2
!      WRITE(10,'(A2,2X,3F20.10)') 'LA',COORDSB(3*(J3-1)+1),COORDSB(3*(J3-1)+2),COORDSB(3*(J3-1)+3)
!  ENDDO
!  CLOSE(10)
!
!  Calculate original centres of mass.
!
CMAX=0.0D0; CMAY=0.0D0; CMAZ=0.0D0
IF (NFREEZE.GT.0) GOTO 11 ! don;t shift or reorient coordinates with frozen atoms
IF (STOCKT.OR.RIGID) THEN 
   DO J1=1,NATOMS/2
      CMAX=CMAX+COORDSA(3*(J1-1)+1)
      CMAY=CMAY+COORDSA(3*(J1-1)+2)
      CMAZ=CMAZ+COORDSA(3*(J1-1)+3)
   ENDDO
   CMAX=2*CMAX/NATOMS; CMAY=2*CMAY/NATOMS; CMAZ=2*CMAZ/NATOMS
   CMBX=0.0D0; CMBY=0.0D0; CMBZ=0.0D0
   DO J1=1,NATOMS/2
      CMBX=CMBX+COORDSB(3*(J1-1)+1)
      CMBY=CMBY+COORDSB(3*(J1-1)+2)
      CMBZ=CMBZ+COORDSB(3*(J1-1)+3)
   ENDDO
   CMBX=2*CMBX/NATOMS; CMBY=2*CMBY/NATOMS; CMBZ=2*CMBZ/NATOMS
ELSE
   DO J1=1,NATOMS
      CMAX=CMAX+COORDSA(3*(J1-1)+1)
      CMAY=CMAY+COORDSA(3*(J1-1)+2)
      CMAZ=CMAZ+COORDSA(3*(J1-1)+3)
   ENDDO
   CMAX=CMAX/NATOMS; CMAY=CMAY/NATOMS; CMAZ=CMAZ/NATOMS
   CMBX=0.0D0; CMBY=0.0D0; CMBZ=0.0D0
   DO J1=1,NATOMS
      CMBX=CMBX+COORDSB(3*(J1-1)+1)
      CMBY=CMBY+COORDSB(3*(J1-1)+2)
      CMBZ=CMBZ+COORDSB(3*(J1-1)+3)
   ENDDO
   CMBX=CMBX/NATOMS; CMBY=CMBY/NATOMS; CMBZ=CMBZ/NATOMS
ENDIF
11 CONTINUE

INVERT=1
DBEST=1.0D100
60 NCHOOSE1=0
65 NCHOOSE1=NCHOOSE1+1
40 NCHOOSE2=0
30 NCHOOSE2=NCHOOSE2+1
DUMMYB(1:3*NATOMS)=COORDSB(1:3*NATOMS)
DUMMYA(1:3*NATOMS)=COORDSA(1:3*NATOMS)
DO J1=1,NATOMS
   ALLPERM(J1)=J1
ENDDO

! The optimal alignment returned by minperdist is a local minimum, but may not
! be the global minimum. Calling MYORIENT first should put permutational isomers
! into a standard alignment and spot the global minimum zedro distance in one
! go. However, we also need to cycle over equivalent atoms in orbits using NCHOOSE2.
!
! Problems can occur if we don't use all the atoms specified by NORBIT1 and NORBIT2
! because of the numerical cutoffs employed in MYORIENT. We could miss the
! right orientation! 
!
! If we use MYORIENT to produce particular orientations then we end up aligning 
! COORDSA not with COORDSB but with the standard orientation of COORDSB in DUMMYB.
! We now deal with this by tracking the complete transformation, including the
! contribution of MYORIENT using ROTB and ROTINVB.
!

IF ((.NOT.BULKT).AND.(NFREEZE.LE.0).AND.(.NOT.CSMT)) THEN
   IF (STOCKT) THEN
      TMAT(1:3,1:3)=0.0D0
      TMAT(1,1)=INVERT*1.0D0; TMAT(2,2)=INVERT*1.0D0; TMAT(3,3)=INVERT*1.0D0
      CALL NEWROTGEOMSTOCK(NATOMS,DUMMYA,TMAT,0.0D0,0.0D0,0.0D0)
      DUMMY(1:3*NATOMS)=DUMMYA(1:3*NATOMS)
      CALL MYORIENT(DUMMYA,DUMMYC,NORBIT1,NCHOOSE1,NORBIT2,NCHOOSE2,NATOMS/2,DEBUG,ROTA,ROTINVA,STOCKT)
      CALL NEWROTGEOMSTOCK(NATOMS,DUMMY,ROTA,0.0D0,0.0D0,0.0D0)
      DUMMYA(1:3*NATOMS)=DUMMY(1:3*NATOMS)

      DUMMY(1:3*NATOMS)=DUMMYB(1:3*NATOMS)
      CALL MYORIENT(DUMMYB,DUMMYC,NORBIT1,1,NORBIT2,1,NATOMS/2,DEBUG,ROTB,ROTINVB,STOCKT)
      CALL NEWROTGEOMSTOCK(NATOMS,DUMMY,ROTB,0.0D0,0.0D0,0.0D0)
      DUMMYB(1:3*NATOMS)=DUMMY(1:3*NATOMS)
      DO J1=1,3*(NATOMS/2)
         DISTANCE=DISTANCE+(DUMMYA(J1)-DUMMYB(J1))**2
      ENDDO
   ELSE
      DUMMYC(1:3*NATOMS)=INVERT*DUMMYA(1:3*NATOMS)
      IF ((PULLT.OR.EFIELDT).AND.(INVERT.EQ.-1)) THEN ! reflect in xz plane
         DO J1=1,NATOMS
            DUMMYC(3*(J1-1)+1)=DUMMYA(3*(J1-1)+1)
            DUMMYC(3*(J1-1)+2)=-DUMMYA(3*(J1-1)+2)
            DUMMYC(3*(J1-1)+3)=DUMMYA(3*(J1-1)+3)
         ENDDO
      ENDIF
      CALL MYORIENT(DUMMYC,DUMMY,NORBIT1,NCHOOSE1,NORBIT2,NCHOOSE2,NATOMS,DEBUG,ROTA,ROTINVA,STOCKT)
      DUMMYA(1:3*NATOMS)=DUMMY(1:3*NATOMS)
      CALL MYORIENT(DUMMYB,DUMMY,NORBIT1,1,NORBIT2,1,NATOMS,DEBUG,ROTB,ROTINVB,STOCKT)
      DUMMYB(1:3*NATOMS)=DUMMY(1:3*NATOMS)
      DISTANCE=0.0D0
      DO J1=1,3*NATOMS
         DISTANCE=DISTANCE+(DUMMYA(J1)-DUMMYB(J1))**2
      ENDDO
   ENDIF
   IF (DEBUG) WRITE(MYUNIT,'(A,G20.10)') ' minpermdist> after initial call to MYORIENT distance=',SQRT(DISTANCE)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  WRITE(10,'(I6)') NATOMS/2
!  WRITE(10,'(A,2I8)') 'DUMMYA after MYORIENT, NCHOOSE1,NCHOOSE2=',NCHOOSE1,NCHOOSE2
!  DO J3=1,NATOMS/2
!      WRITE(10,'(A2,2X,3F20.10)') 'LA',DUMMYA(3*(J3-1)+1),DUMMYA(3*(J3-1)+2),DUMMYA(3*(J3-1)+3)
!  ENDDO
!  WRITE(10,'(I6)') NATOMS/2
!  WRITE(10,'(A)') 'DUMMYB after MYORIENT'
!  DO J3=1,NATOMS/2
!      WRITE(10,'(A2,2X,3F20.10)') 'LA',DUMMYB(3*(J3-1)+1),DUMMYB(3*(J3-1)+2),DUMMYB(3*(J3-1)+3)
!  ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ELSE
  NORBIT1=1
  NORBIT2=1
  ROTINVB(1:3,1:3)=0.0D0
  ROTINVB(1,1)=1.0D0; ROTINVB(2,2)=1.0D0; ROTINVB(3,3)=1.0D0
  CALL NEWMINDIST(DUMMYB,DUMMYA,NATOMS,DISTANCE,BULKT,TWOD,'AX    ',.FALSE.,RIGID,DEBUG,RMAT)
  IF (DEBUG) WRITE(MYUNIT,'(A,G20.10)') ' minpermdist> after initial call to NEWMINDIST distance=',DISTANCE
  DISTANCE=DISTANCE**2 ! minpermdist returns the distance squared for historical reasons
ENDIF
IF (DEBUG) WRITE(MYUNIT,'(A,4I8)') ' minpermdist> size of orbits and selected atoms: ',NORBIT1,NORBIT2,NCHOOSE1,NCHOOSE2

! CALL POTENTIAL(DUMMYA,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
! PRINT '(2(A,F25.15))',' New initial energy=',ENERGY,' RMS=',RMS
! PRINT '(2(A,F25.15))',' for coordinates:'
! PRINT '(3F25.15)',DUMMYA(1:3*NATOMS)
!
!  Bipartite matching routine for permutations. Coordinates in DUMMYB do not change
!  but the coordinates in DUMMYA do. DISTANCE is the distance in this case.
!  We return to label 10 after every round of permutational/orientational alignment
!  unless we have converged to the identity permutation.
!
!  Atoms are not allowed to appear in more than one group.
!  The maximum number of pair exchanges associated with a group is two.
!
NTRIES=0
!
!  RMATCUMUL contains the accumulated rotation matrix that relates the original 
!  DUMMYA obtained from COORDSA to the final one.
!
RMATCUMUL(1:3,1:3)=0.0D0
RMATCUMUL(1,1)=1.0D0; RMATCUMUL(2,2)=1.0D0; RMATCUMUL(3,3)=1.0D0
PDISTANCE=1.0D100
10 CONTINUE
NTRIES=NTRIES+1

NDUMMY=1
DO J1=1,NATOMS
   NEWPERM(J1)=J1
ENDDO
!
! ALLPERM saves the permutation from the previous cycle.
! NEWPERM contains the permutation for this cycle, relative to the identity.
! SAVEPERM is temporary storage for NEWPERM.
! NEWPERM must be applied to ALLPERM after the loop over NPERMGROUP and
! corresponding swaps.
!
! New version allows for overlapping atoms in NPERMGROUP, so that atoms
! can appear in more thsan one group. This was needed to flexible water potentials.
!
DO J1=1,NPERMGROUP
   PATOMS=NPERMSIZE(J1)
   DO J2=1,PATOMS
      PDUMMYA(3*(J2-1)+1)=DUMMYA(3*(NEWPERM(PERMGROUP(NDUMMY+J2-1))-1)+1)
      PDUMMYA(3*(J2-1)+2)=DUMMYA(3*(NEWPERM(PERMGROUP(NDUMMY+J2-1))-1)+2)
      PDUMMYA(3*(J2-1)+3)=DUMMYA(3*(NEWPERM(PERMGROUP(NDUMMY+J2-1))-1)+3)
      PDUMMYB(3*(J2-1)+1)=DUMMYB(3*(NEWPERM(PERMGROUP(NDUMMY+J2-1))-1)+1)
      PDUMMYB(3*(J2-1)+2)=DUMMYB(3*(NEWPERM(PERMGROUP(NDUMMY+J2-1))-1)+2)
      PDUMMYB(3*(J2-1)+3)=DUMMYB(3*(NEWPERM(PERMGROUP(NDUMMY+J2-1))-1)+3)
   ENDDO
   CALL MINPERM(PATOMS, PDUMMYB, PDUMMYA, BOXLX, BOXLY, BOXLZ, BULKT, LPERM, LDISTANCE, DIST2, WORSTRAD)
   SAVEPERM(1:NATOMS)=NEWPERM(1:NATOMS)
   DO J2=1,PATOMS
      SAVEPERM(PERMGROUP(NDUMMY+J2-1))=NEWPERM(PERMGROUP(NDUMMY+LPERM(J2)-1))
   ENDDO
!
! Update permutation of associated atoms, if any.
! We must do this as we go along, because these atoms could move in more than
! one permutational group now.
!
   IF (NSETS(J1).GT.0) THEN
      DO J2=1,PATOMS
         DO J3=1,NSETS(J1)
            SAVEPERM(SETS(PERMGROUP(NDUMMY+J2-1),J3))=SETS(NEWPERM(PERMGROUP(NDUMMY+LPERM(J2)-1)),J3)
         ENDDO
      ENDDO
   ENDIF
   NDUMMY=NDUMMY+NPERMSIZE(J1)
   NEWPERM(1:NATOMS)=SAVEPERM(1:NATOMS)
ENDDO
!
! Update the overall permutation here.
! 
DO J1=1,NATOMS
   SAVEPERM(ALLPERM(J1))=ALLPERM(NEWPERM(J1))
ENDDO
ALLPERM(1:NATOMS)=SAVEPERM(1:NATOMS)

DUMMY(1:3*NATOMS)=DUMMYA(1:3*NATOMS)
NPERM=0
DISTANCE=0.0D0
IF (STOCKT .OR. RIGID) THEN ! additional permutation of dipoles
   DO J1=(NATOMS/2)+1,NATOMS
      ALLPERM(J1)=ALLPERM(J1-(NATOMS/2))+(NATOMS/2)
      NEWPERM(J1)=NEWPERM(J1-(NATOMS/2))+(NATOMS/2)
   ENDDO
ENDIF
!
! Update coordinates in DUMMYA to overall permutation using NEWPERM.
!
DO J3=1,NATOMS
   DUMMYA(3*(J3-1)+1)=DUMMY(3*(NEWPERM(J3)-1)+1)
   DUMMYA(3*(J3-1)+2)=DUMMY(3*(NEWPERM(J3)-1)+2)
   DUMMYA(3*(J3-1)+3)=DUMMY(3*(NEWPERM(J3)-1)+3)
   IF (J3.NE.NEWPERM(J3)) THEN
!     IF (DEBUG) WRITE(*,'(A,I5,A,I5)') ' minpermdist> move position ',NEWPERM(J3),' to ',J3
      NPERM=NPERM+1
   ENDIF
   IF (STOCKT .OR. RIGID) THEN
      IF (J3.LE.(NATOMS/2)) THEN
         DISTANCE=DISTANCE+(DUMMYA(3*(J3-1)+1)-DUMMYB(3*(J3-1)+1))**2 &
  &                       +(DUMMYA(3*(J3-1)+2)-DUMMYB(3*(J3-1)+2))**2 &
  &                       +(DUMMYA(3*(J3-1)+3)-DUMMYB(3*(J3-1)+3))**2
      ENDIF
   ELSEIF (.NOT.BULKT) THEN
      DISTANCE=DISTANCE+(DUMMYA(3*(J3-1)+1)-DUMMYB(3*(J3-1)+1))**2 &
  &                    +(DUMMYA(3*(J3-1)+2)-DUMMYB(3*(J3-1)+2))**2 &
  &                    +(DUMMYA(3*(J3-1)+3)-DUMMYB(3*(J3-1)+3))**2
   ELSE
      DISTANCE=DISTANCE + (DUMMYA(3*(J3-1)+1)-DUMMYB(3*(J3-1)+1)- BOXLX*NINT((DUMMYA(3*(J3-1)+1)-DUMMYB(3*(J3-1)+1))/BOXLX))**2 &
  &                     + (DUMMYA(3*(J3-1)+2)-DUMMYB(3*(J3-1)+2)- BOXLY*NINT((DUMMYA(3*(J3-1)+2)-DUMMYB(3*(J3-1)+2))/BOXLY))**2
      IF (.NOT.TWOD) DISTANCE=DISTANCE+(DUMMYA(3*(J3-1)+3)-DUMMYB(3*(J3-1)+3) -  &
  &                                                               BOXLZ*NINT((DUMMYA(3*(J3-1)+3)-DUMMYB(3*(J3-1)+3))/BOXLZ))**2
   ENDIF
ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALL UPDATENBONDS(DUMMYA)
! CALL OCHARMM(DUMMYA,VNEW,ENERGY,.FALSE.,.FALSE.)
! CALL POTENTIAL(DUMMYA,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
! PRINT '(2(A,F25.15))',' B Energy is now=',ENERGY,' RMS=',RMS
! PRINT '(2(A,F25.15))',' for coordinates:'
! PRINT '(3F25.15)',DUMMYA(1:3*NATOMS)
! PRINT '(A,F25.15,A)',' Energy is now=',ENERGY,' kcal/mol'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF (DEBUG) WRITE(MYUNIT,'(A,I6,A,G20.10)') ' minpermdist> distance after permuting ',NPERM,' pairs of atoms=',SQRT(DISTANCE)

! CALL OCHARMM(DUMMYA,VNEW,ENERGY,.FALSE.,.FALSE.)
! PRINT '(A,F25.15,A)',' Energy for last cycle=',ENERGY,' kcal/mol'
! CALL UPDATENBONDS(DUMMYA)
! PRINT '(A,F25.15,A)',' Energy for last cycle=',ENERGY,' kcal/mol after update'
! CALL POTENTIAL(DUMMYA,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
! PRINT '(2(A,F25.15))',' Energy for last cycle=',ENERGY,' RMS=',RMS
!
!  Optimal alignment. Coordinates in DUMMYA are reset by NEWMINDIST (second argument).
!  Must allow at least one call to NEWMINDIST in case the MYORIENT result is terrible
!  but gives zero permutations!
!  
!          IF (DEBUG) THEN
!             OPEN(UNIT=765,FILE='DUMMYB.xyz',STATUS='UNKNOWN')
!             NDUMMY=(NATOMS/CSMGPINDEX)
!             DO J1=1,CSMGPINDEX
!                WRITE(765,'(I6)') NDUMMY
!                WRITE(765,'(A,I6)') 'result of point group operation ',J1
!                WRITE(765,'(A,3G20.10)') ('LA ',DUMMYB(3*(J2-1)+1+3*NDUMMY*(J1-1)), &
!      &                    DUMMYB(3*(J2-1)+2+3*NDUMMY*(J1-1)),DUMMYB(3*(J2-1)+3+3*NDUMMY*(J1-1)),J2=1,NDUMMY)
!             ENDDO
!             CLOSE(765)
!             OPEN(UNIT=765,FILE='DUMMYA.xyz',STATUS='UNKNOWN')
!             DO J1=1,CSMGPINDEX
!                WRITE(765,'(I6)') NDUMMY
!                WRITE(765,'(A,I6)') 'matching point group operation ',J1
!                WRITE(765,'(A,3G20.10)') ('LA ',DUMMYA(3*(J2-1)+1+3*NDUMMY*(J1-1)), &
!      &                    DUMMYA(3*(J2-1)+2+3*NDUMMY*(J1-1)),DUMMYA(3*(J2-1)+3+3*NDUMMY*(J1-1)),J2=1,NDUMMY)
!             ENDDO
!             CLOSE(765)
!          ENDIF
! STOP

IF ((NPERM.NE.0).OR.(NTRIES.EQ.1)) THEN 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  PRINT '(A)','DUMMYA before NEWMINDIST:'
!  PRINT '(3F20.10)',DUMMYA(1:3*(NATOMS/2))
!  PRINT '(A)','DUMMYB before NEWMINDIST:'
!  PRINT '(3F20.10)',DUMMYB(1:3*(NATOMS/2))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   CALL NEWMINDIST(DUMMYB,DUMMYA,NATOMS,DISTANCE,BULKT,TWOD,'AX    ',.FALSE.,RIGID,DEBUG,RMAT)
   RMATCUMUL=MATMUL(RMAT,RMATCUMUL)
   DISTANCE=DISTANCE**2 ! we are using DISTANCE^2 further down

   IF (DEBUG) WRITE(MYUNIT,'(A,G20.10,A,I6)')  ' minpermdist> distance after NEWMINDIST=                     ',SQRT(DISTANCE), &
  &           ' tries=',NTRIES  
   CALL FLUSH(MYUNIT,ISTAT)
   IF ((NTRIES.LT.MAXIMUMTRIES).AND.(DISTANCE.LT.PDISTANCE)) THEN
!     WRITE(MYUNIT,'(A,I6,2G20.10)') 'NTRIES,DISTANCE,PDISTANCE=',NTRIES,DISTANCE,PDISTANCE
      PDISTANCE=DISTANCE
      GOTO 10
   ELSE ! prevent infinite loop
      IF (DEBUG) WRITE(MYUNIT,'(A)') &
  &              ' minpermdist> WARNING - number of tries exceeded or distance increased with nonzero permutations'
   ENDIF
ENDIF

IF (DISTANCE.LT.DBEST) THEN
   DBEST=DISTANCE
   XBEST(1:3*NATOMS)=DUMMYA(1:3*NATOMS)
   RMATBEST(1:3,1:3)=RMATCUMUL(1:3,1:3)
   ROTINVBBEST(1:3,1:3)=ROTINVB(1:3,1:3) 
   ROTABEST(1:3,1:3)=ROTA(1:3,1:3)      
   RMATBEST=MATMUL(RMATBEST,ROTABEST)
   IF (INVERT.EQ.-1) THEN
      IF (PULLT.OR.EFIELDT) THEN ! reflect in xz plane rather than invert!
         RMATBEST(1:3,1:3)=MATMUL(RMATBEST,REFXZ)
      ELSE
         RMATBEST(1:3,1:3)=-RMATBEST(1:3,1:3)
      ENDIF
   ENDIF
ENDIF

!
! If GEOMDIFFTOL is set too small we could miss the best solution by exiting prematurely. 
! Turn off the next line?!
!
! IF (SQRT(DBEST).LT.GEOMDIFFTOL) GOTO 50
IF (NCHOOSE2.LT.NORBIT2) GOTO 30
IF (NCHOOSE1.LT.NORBIT1) GOTO 65
!
!  Now try the enantiomer (or xz reflected structure for PULLT.OR.EFIELDT).
!
IF ((NCHOOSE2.EQ.NORBIT2).AND.(NCHOOSE1.EQ.NORBIT1).AND.(INVERT.EQ.1)) THEN
!
! don't try inversion for bulk or charmm or amber or frozen atoms or CSM
!
   IF (BULKT.OR.CHRMMT.OR.AMBERT.OR.(NFREEZE.GT.0).OR.CSMT) GOTO 50 
   IF (DEBUG) WRITE(MYUNIT,'(A)') ' minpermdist> inverting geometry for comparison with target'
   INVERT=-1
   GOTO 60
ENDIF

50 DISTANCE=DBEST
!
!  XBEST contains the best alignment of A coordinates for the orientation of B coordinates in DUMMYB.
!  Rotate XBEST by ROTINVB to put in best correspondence with COORDSB, undoing the reorientation to DUMMYB from MYORIENT. 
!  We should get the same result for ROTINVB * RMATBEST * (COORDSA-CMA) 
!  where RMATBEST = +/- RMATCUMUL * ROTA for the best alignment 
!  (aside from a possible permutation of the atom ordering)
!
   IF (NFREEZE.GT.0) THEN
      XDUMMY=0.0D0
      DO J1=1,NATOMS
         XDUMMY=XDUMMY+(COORDSB(3*(J1-1)+1)-XBEST(3*(J1-1)+1))**2+ &
  &                    (COORDSB(3*(J1-1)+2)-XBEST(3*(J1-1)+2))**2+ &
  &                    (COORDSB(3*(J1-1)+3)-XBEST(3*(J1-1)+3))**2
      ENDDO
   ELSE IF (STOCKT) THEN
      CALL NEWROTGEOMSTOCK(NATOMS,XBEST,ROTINVBBEST,0.0D0,0.0D0,0.0D0)
      XDUMMY=0.0D0
      DO J1=1,(NATOMS/2)
         XBEST(3*(J1-1)+1)=XBEST(3*(J1-1)+1)+CMBX
         XBEST(3*(J1-1)+2)=XBEST(3*(J1-1)+2)+CMBY
         XBEST(3*(J1-1)+3)=XBEST(3*(J1-1)+3)+CMBZ
         XDUMMY=XDUMMY+(COORDSB(3*(J1-1)+1)-XBEST(3*(J1-1)+1))**2+ &
  &                    (COORDSB(3*(J1-1)+2)-XBEST(3*(J1-1)+2))**2+ &
  &                    (COORDSB(3*(J1-1)+3)-XBEST(3*(J1-1)+3))**2
      ENDDO
   ELSEIF (BULKT) THEN
      XDUMMY=0.0D0
      DO J1=1,NATOMS
         XDUMMY=XDUMMY+(COORDSB(3*(J1-1)+1)-XBEST(3*(J1-1)+1) - BOXLX*NINT((COORDSB(3*(J1-1)+1)-XBEST(3*(J1-1)+1))/BOXLX))**2+ &
  &                    (COORDSB(3*(J1-1)+2)-XBEST(3*(J1-1)+2) - BOXLY*NINT((COORDSB(3*(J1-1)+2)-XBEST(3*(J1-1)+2))/BOXLY))**2
         IF (.NOT.TWOD) XDUMMY=XDUMMY+(COORDSB(3*(J1-1)+3)-XBEST(3*(J1-1)+3) - &
  &                                                             BOXLZ*NINT((COORDSB(3*(J1-1)+3)-XBEST(3*(J1-1)+3))/BOXLZ))**2
      ENDDO
   ELSE
      XDUMMY=0.0D0
      DO J1=1,NATOMS
         XBEST(3*(J1-1)+1:3*(J1-1)+3)=MATMUL(ROTINVBBEST,XBEST(3*(J1-1)+1:3*(J1-1)+3))
         XBEST(3*(J1-1)+1)=XBEST(3*(J1-1)+1)+CMBX
         XBEST(3*(J1-1)+2)=XBEST(3*(J1-1)+2)+CMBY
         XBEST(3*(J1-1)+3)=XBEST(3*(J1-1)+3)+CMBZ
         XDUMMY=XDUMMY+(COORDSB(3*(J1-1)+1)-XBEST(3*(J1-1)+1))**2+ &
  &                    (COORDSB(3*(J1-1)+2)-XBEST(3*(J1-1)+2))**2+ &
  &                    (COORDSB(3*(J1-1)+3)-XBEST(3*(J1-1)+3))**2
      ENDDO
   ENDIF
   IF (ABS(SQRT(XDUMMY)-SQRT(DISTANCE)).GT.GEOMDIFFTOL) THEN
      WRITE(MYUNIT,'(2(A,F20.10))') ' minpermdist> ERROR *** distance between transformed XBEST and COORDSB=',SQRT(XDUMMY), &
  &                         ' should be ',SQRT(DISTANCE)
   ENDIF

   IF (NFREEZE.GT.0) THEN
      RMATBEST(1:3,1:3)=0.0D0
      RMATBEST(1,1)=1.0D0; RMATBEST(2,2)=1.0D0; RMATBEST(3,3)=1.0D0
   ELSE
      RMATBEST=MATMUL(ROTINVB,RMATBEST)
   ENDIF
   COORDSA(1:3*NATOMS)=XBEST(1:3*NATOMS) ! finally, best COORDSA should include permutations for DNEB input!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  DO J1=1,(NATOMS/2)
!     XDUMMY=XDUMMY+(COORDSB(3*(J1-1)+1)-COORDSA(3*(J1-1)+1))**2+ &
! &                 (COORDSB(3*(J1-1)+2)-COORDSA(3*(J1-1)+2))**2+ &
! &                 (COORDSB(3*(J1-1)+3)-COORDSA(3*(J1-1)+3))**2
!  ENDDO
!  PRINT '(A,F20.10)','XDUMMY=',XDUMMY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DISTANCE=SQRT(DISTANCE)

RETURN
END SUBROUTINE MINPERMDIST
