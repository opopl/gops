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

USE KEY,ONLY : NPERMGROUP, NPERMSIZE, PERMGROUP, NSETS, SETS, STOCKT, GEOMDIFFTOL, AMBERT, &
  &            NFREEZE, NABT, RBAAT, ANGLEAXIS2, BESTPERM, LOCALPERMDIST, PULLT, EFIELDT, NTSITES, RIGIDBODY, PERMDIST
USE MODCHARMM,ONLY : CHRMMT
USE MODAMBER9, ONLY: NOPERMPROCHIRAL, PROCHIRALH
USE INTCOMMONS, ONLY : INTMINPERMT, INTINTERPT, DESMINT, OLDINTMINPERMT, INTDISTANCET
USE INTCUTILS, ONLY : INTMINPERM, OLD_INTMINPERM, INTMINPERM_CHIRAL, INTDISTANCE
IMPLICIT NONE

INTEGER, PARAMETER :: MAXIMUMTRIES=100
INTEGER NATOMS, NPERM, PATOMS, NTRIES, NRB
INTEGER J3, J4, INVERT, NORBIT1, NORBIT2, NCHOOSE2, NDUMMY, LPERM(NATOMS), J1, J2, NCHOOSE1
DOUBLE PRECISION DIST2, COORDSA(3*NATOMS), COORDSB(3*NATOMS), DISTANCE, DUMMYA(3*NATOMS), &
  &              BESTA(3*NATOMS), DUMMYB(3*NATOMS), DUMMY(3*NATOMS)
DOUBLE PRECISION BOXLX,BOXLY,BOXLZ,WORSTRAD,RMAT(3,3),ENERGY, VNEW(3*NATOMS), DX, DY, DZ, RMS, DBEST, XBEST(3*NATOMS)
DOUBLE PRECISION CMXA, CMXB, CMXC, QBEST(4), SITESA(3*NTSITES), SITESB(3*NTSITES)
DOUBLE PRECISION ROTA(3,3), ROTINVA(3,3), ROTB(3,3), ROTINVB(3,3), ROTINVBBEST(3,3), ROTABEST(3,3), RMATBEST(3,3), TMAT(3,3)
DOUBLE PRECISION CMAX, CMAY, CMAZ, CMBX, CMBY, CMBZ, RMATCUMUL(3,3), PVEC(3), RTEMP1(3,3), RTEMP2(3,3)
DOUBLE PRECISION REFXZ(3,3)
LOGICAL DEBUG, TWOD, RIGID, BULKT, PITEST
DOUBLE PRECISION PDUMMYA(3*NATOMS), PDUMMYB(3*NATOMS), LDISTANCE, DUMMYC(3*NATOMS), XDUMMY, DUMMYD(3*NATOMS), LDBEST
INTEGER NEWPERM(NATOMS), ALLPERM(NATOMS), SAVEPERM(NATOMS)
DOUBLE PRECISION TIME0, TIME1, BSX, BSY, BSZ
DOUBLE PRECISION, ALLOCATABLE :: TEMPA(:), TEMPB(:)
COMMON /BULKSHIFT/ BSX, BSY, BSZ

! CALL POTENTIAL(COORDSA,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
! PRINT '(2(A,F25.15))',' intial energy for A=             ',ENERGY,' RMS=',RMS
! CALL POTENTIAL(COORDSB,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
! PRINT '(2(A,F25.15))',' intial energy for B=             ',ENERGY,' RMS=',RMS
!
! For angle-axis coordinates with PERMDIST: 
! (1) use MINPERM to permute the centre-of-mass coordinates in the usual way,
!     using a metric based on just the centre of mass. These coordinates are
!     stored in the first 3*NATOMS entries.
! (2) for each reorientation of the centre of mass corrdinates we have to
!     rotate the orientational coordinates. Then we need a loop over the
!     rigid bodies to minimise the distance metric based upon all the sites
!     for the allowed internal symmetry operations of every rigid body.
!

  IF (RBAAT) THEN
     CALL RBMINPERMDIST(COORDSB,COORDSA,DISTANCE,DIST2,QBEST,RMATBEST,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,SITESB,SITESA)
     RETURN
  ENDIF
!
REFXZ(1:3,1:3)=0.0D0
REFXZ(1,1)=1.0D0; REFXZ(2,2)=-1.0D0; REFXZ(3,3)=1.0D0
!
! The INTMINPERM keyword may now be OK for movie making. It was producing jumps
! but I think setting RMATBEST has fixed this. DJW 18/2/11.
!

IF (NOPERMPROCHIRAL.AND.(AMBERT.OR.NABT).AND..NOT.INTMINPERMT) THEN
   CALL NEWMINDIST(COORDSB,COORDSA,NATOMS,DISTANCE,BULKT,TWOD,'AX   ',.FALSE.,RIGID,DEBUG,RMAT)
   IF (.NOT.ALLOCATED(PROCHIRALH)) CALL FINDCHIRALH(DEBUG)
   CALL MINPERM_CHIRAL(COORDSB, COORDSA,DISTANCE,RMAT, DEBUG)
   CALL check_valleu_chirality(COORDSB, COORDSA,DEBUG)
   CALL NEWMINDIST(COORDSB,COORDSA,NATOMS,DISTANCE,BULKT,TWOD,'AX   ',.TRUE.,RIGID,DEBUG,RMAT)
   RMATBEST = RMAT 
   RETURN
ENDIF

IF (.NOT.PERMDIST) THEN
!
! NEWMINDIST is always called with PRESERVET .FALSE. here. Hence COORDSA will generally be
! changed to be put into best correspondence with COORDSB.
!
   CALL NEWMINDIST(COORDSB,COORDSA,NATOMS,DISTANCE,BULKT,TWOD,'AX   ',.FALSE.,RIGIDBODY,DEBUG,RMAT)
   RETURN
ENDIF

IF (INTMINPERMT.AND.(INTINTERPT.OR.DESMINT)) THEN
    IF (CHRMMT.OR.OLDINTMINPERMT) THEN
      !CALL MYCPU_TIME(TIME0,.FALSE.)
      CALL OLD_INTMINPERM(COORDSB, COORDSA, DISTANCE, RMAT, DEBUG)
      CALL NEWMINDIST(COORDSB,COORDSA,NATOMS,DISTANCE,BULKT,TWOD,'AX    ',.FALSE.,RIGID,DEBUG,RMAT)
      !CALL MYCPU_TIME(TIME1,.FALSE.)
      !PRINT*, "alignment took", TIME1-TIME0
      IF (DEBUG) PRINT '(A,G20.10)', "minpermdist> newmin-distance ", DISTANCE
      DISTANCE = DISTANCE**2
      RMATBEST=RMAT
    ELSEIF (AMBERT.OR.NABT) THEN
      IF (NOPERMPROCHIRAL) THEN
        IF (.NOT.ALLOCATED(PROCHIRALH)) CALL FINDCHIRALH(DEBUG)
        CALL INTMINPERM_CHIRAL(COORDSB, COORDSA, DISTANCE, RMAT, DEBUG)
      ELSE
        CALL INTMINPERM(COORDSB,COORDSA,DISTANCE,RMAT,DEBUG)
      ENDIF
      IF (DEBUG) PRINT*, "distance in minperm", SQRT(DISTANCE)
      CALL check_valleu_chirality(COORDSB, COORDSA,DEBUG)
      CALL NEWMINDIST(COORDSB,COORDSA,NATOMS,DISTANCE,BULKT,TWOD,'AX   ',.FALSE.,RIGID,DEBUG,RMAT)
      IF (DEBUG) PRINT '(A,G20.10)',"minpermdist> newmin-distance ", DISTANCE
      DISTANCE = DISTANCE**2 ! see below
      RMATBEST=RMAT
    ELSE
      PRINT*, "minpermdist> ERROR *** using INTMINPERM without CHARMM/AMBER"
      STOP
    ENDIF
    IF (INTDISTANCET) THEN
      CALL INTDISTANCE(COORDSB, COORDSA, DISTANCE, DEBUG)
      IF (DEBUG) PRINT*, "msb50 minpermdist using intdistance", DISTANCE
    ENDIF
   DISTANCE=SQRT(DISTANCE)
   RETURN
ENDIF
!
!  Calculate original centres of mass.
!
CMAX=0.0D0; CMAY=0.0D0; CMAZ=0.0D0
IF (NFREEZE.GT.0) GOTO 11 ! don;t shift or reorient coordinates with frozen atoms
IF (STOCKT) THEN 
   NRB=(NATOMS/2)
   DO J1=1,NRB
      CMAX=CMAX+COORDSA(3*(J1-1)+1)
      CMAY=CMAY+COORDSA(3*(J1-1)+2)
      CMAZ=CMAZ+COORDSA(3*(J1-1)+3)
   ENDDO
   CMAX=CMAX/NRB; CMAY=CMAY/NRB; CMAZ=CMAZ/NRB
   CMBX=0.0D0; CMBY=0.0D0; CMBZ=0.0D0
   DO J1=1,NRB
      CMBX=CMBX+COORDSB(3*(J1-1)+1)
      CMBY=CMBY+COORDSB(3*(J1-1)+2)
      CMBZ=CMBZ+COORDSB(3*(J1-1)+3)
   ENDDO
   CMBX=CMBX/NRB; CMBY=CMBY/NRB; CMBZ=CMBZ/NRB
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

DISTANCE=0.0D0
IF (NFREEZE.LE.0) THEN
   IF (BULKT) THEN
      NORBIT1=1; NORBIT2=1
      CALL BULKMINDIST(DUMMYB,DUMMYA,NATOMS,DISTANCE,TWOD,DEBUG,BOXLX,BOXLY,BOXLZ,PITEST,.TRUE.)
      IF (PITEST) THEN
         COORDSA(1:3*NATOMS)=DUMMYA(1:3*NATOMS)
         RMATBEST(1:3,1:3)=0.0D0
         RMATBEST(1,1)=1.0D0; RMATBEST(2,2)=1.0D0; RMATBEST(3,3)=1.0D0
         DISTANCE=SQRT(DISTANCE)
         RETURN
      ELSE
         CALL NEWMINDIST(DUMMYB,DUMMYA,NATOMS,DISTANCE,BULKT,TWOD,'AX    ',.FALSE.,RIGID,DEBUG,RMAT)
         IF (DEBUG) PRINT '(A,G20.10)',' minpermdist> after initial call to BULK/NEWMINDIST distance=',DISTANCE
         DISTANCE=DISTANCE**2 ! minperdist returns the distance squared for historical reasons
      ENDIF
   ELSEIF (STOCKT) THEN
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

!     PRINT *, 'DUMMYA'
!     PRINT *, DUMMYA
!     PRINT *, 'DUMMYB'
!     PRINT *, DUMMYB
!      STOP

!  IF (DEBUG) PRINT '(A,G20.10,A,I6,A)', &
! &       ' minpermdist> after initial call to MYORIENT distance=',SQRT(DISTANCE), ' for ',NATOMS,' atoms'
!  IF (DEBUG) PRINT '(A,4I8)',' minpermdist> size of orbits and selected atoms: ', &
! &       NORBIT1,NORBIT2,NCHOOSE1,NCHOOSE2
ELSE
   NORBIT1=1; NORBIT2=1
   CALL NEWMINDIST(DUMMYB,DUMMYA,NATOMS,DISTANCE,BULKT,TWOD,'AX    ',.FALSE.,RIGID,DEBUG,RMAT)
   IF (DEBUG) PRINT '(A,G20.10)',' minpermdist> after initial call to NEWMINDIST distance=',DISTANCE
   DISTANCE=DISTANCE**2
ENDIF
!
!  Bipartite matching routine for permutations. Coordinates in DUMMYB do not change
!  but the coordinates in DUMMYA do. DISTANCE is the distance^2 in this case.
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
!
! All permutations within this group of size NPERMSIZE(J1) are now tried.
!
   CALL MINPERM(PATOMS, PDUMMYB, PDUMMYA, BOXLX, BOXLY, BOXLZ, BULKT, LPERM, LDISTANCE, DIST2, WORSTRAD)
   SAVEPERM(1:NATOMS)=NEWPERM(1:NATOMS)
   DO J2=1,PATOMS
      saveperm(permgroup(ndummy+j2-1))=newperm(permgroup(ndummy+lperm(j2)-1))
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
IF (STOCKT) THEN ! additional permutation of dipoles.
   DO J1=(NATOMS/2)+1,NATOMS
      ALLPERM(J1)=ALLPERM(J1-(NATOMS/2))+(NATOMS/2)
      NEWPERM(J1)=NEWPERM(J1-(NATOMS/2))+(NATOMS/2)
   ENDDO
ELSE IF (ANGLEAXIS2) THEN ! additional permutation for angle-axis variables - this is the obsolete ANGLEAXIS!
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
   IF (STOCKT.OR.ANGLEAXIS2) THEN
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

! CALL POTENTIAL(DUMMYA,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
! PRINT '(2(A,F25.15))',' energy for A now          ',ENERGY,' RMS=',RMS

! IF (DEBUG) WRITE(*,'(A,I6,A,G20.10)') ' minpermdist> distance after moving ',NPERM,' atoms=',SQRT(DISTANCE)

! CALL OCHARMM(DUMMYA,VNEW,ENERGY,.FALSE.,.FALSE.)
! PRINT '(A,F25.15,A)',' Energy=',ENERGY,' kcal/mol'
! CALL UPDATENBONDS(DUMMYA)
! PRINT '(A,F25.15,A)',' Energy=',ENERGY,' kcal/mol after update'
! WRITE(*,'(A,I6,A,G20.10)') ' minpermdist> distance after permuting ',NPERM,' pairs of atoms=',SQRT(DISTANCE)
!
!  Optimal alignment. Coordinates in DUMMYA are reset by NEWMINDIST (second argument).
!  Must allow at least one call to NEWMINDIST in case the MYORIENT result is terrible
!  but gives zero permutations!
!  
 
! PRINT '(A,I6,2G20.10)','NPERM,DBEST,DISTANCE=',NPERM,DBEST,DISTANCE
IF ((NPERM.NE.0).OR.(NTRIES.EQ.1)) THEN 
   CALL NEWMINDIST(DUMMYB,DUMMYA,NATOMS,DISTANCE,BULKT,TWOD,'AX    ',.FALSE.,RIGID,DEBUG,RMAT)
   RMATCUMUL=MATMUL(RMAT,RMATCUMUL)
   DISTANCE=DISTANCE**2 ! we are using DISTANCE^2 further down
!  IF (DEBUG) WRITE(*,'(A,G20.10)') ' minpermdist> distance after NEWMINDIST=                     ', &
! &                                    SQRT(DISTANCE) 
   IF (NTRIES.LT.MAXIMUMTRIES) THEN
      GOTO 10
   ELSE ! prevent infinite loop
      PRINT '(A)',' minpermdist> WARNING - number of tries exceeded, giving up'
   ENDIF
ENDIF

IF (DISTANCE.LT.DBEST) THEN
   DBEST=DISTANCE
!  PRINT *, 'DBEST=', DBEST
   XBEST(1:3*NATOMS)=DUMMYA(1:3*NATOMS)
!  PRINT *, 'XBEST'
!  PRINT *, XBEST
   BESTPERM(1:NATOMS)=ALLPERM(1:NATOMS)
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
!  GO TO 50
IF ((NCHOOSE2.EQ.NORBIT2).AND.(NCHOOSE1.EQ.NORBIT1).AND.(INVERT.EQ.1)) THEN
!
! don't try inversion for bulk or charmm or amber or frozen atoms
!
   IF (BULKT.OR.CHRMMT.OR.AMBERT.OR.NABT.OR.(NFREEZE.GT.0)) GOTO 50 
  IF (DEBUG) PRINT '(A)',' minpermdist> inverting geometry for comparison with target'
   INVERT=-1
   GOTO 60
ENDIF

50 DISTANCE=DBEST
!
!  XBEST contains the best alignment of A coordinates for the orientation of B coordinates in DUMMYB.
!  Rotate XBEST by ROTINVBBEST to put in best correspondence with COORDSB, undoing the reorientation to DUMMYB from MYORIENT. 
!  We should get the same result for ROTINVBBEST * RMATBEST * (COORDSA-CMA) 
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
!     PRINT *, CMBX, CMBY, CMBZ
!     PRINT *, 'ROTINVBBEST'
!     PRINT *, ROTINVBBEST
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
      PRINT '(2(A,F20.10))',' minpermdist> ERROR *** distance between transformed XBEST and COORDSB=',SQRT(XDUMMY), &
  &                         ' should be ',SQRT(DISTANCE)
      STOP
   ENDIF

   IF (NFREEZE.GT.0) THEN
      RMATBEST(1:3,1:3)=0.0D0
      RMATBEST(1,1)=1.0D0; RMATBEST(2,2)=1.0D0; RMATBEST(3,3)=1.0D0
   ELSE
      RMATBEST=MATMUL(ROTINVBBEST,RMATBEST)
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
!     CLOSE(10)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALL POTENTIAL(COORDSA,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
! PRINT '(2(A,F25.15))',' before check_valleau A=',ENERGY,' RMS=',RMS
! CALL POTENTIAL(COORDSB,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
! PRINT '(2(A,F25.15))',' before check_valleau B=',ENERGY,' RMS=',RMS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      IF ((AMBERT.OR.NABT).AND.(.NOT.LOCALPERMDIST)) THEN
         CALL check_valleu_chirality(COORDSB, COORDSA,DEBUG)
         CALL NEWMINDIST(COORDSB,COORDSA,NATOMS,DISTANCE,BULKT,TWOD,'AX    ',.TRUE.,RIGID,DEBUG,RMAT)
         DISTANCE=DISTANCE**2 ! minpermdist used to return the distance squared for historical reasons!
      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALL POTENTIAL(COORDSA,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
! PRINT '(2(A,F25.15))',' before check_valleau A=',ENERGY,' RMS=',RMS
! CALL POTENTIAL(COORDSB,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
! PRINT '(2(A,F25.15))',' before check_valleau B=',ENERGY,' RMS=',RMS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! IF (DEBUG) PRINT '(A)',' minpermdist> Overall permutation for COORDSA (second argument):'
! IF (DEBUG) PRINT '(20I6)',BESTPERM(1:NATOMS)
! PRINT '(I6)',NATOMS
! PRINT '(A)','coordsa in minpermdist:'
! PRINT '(A,3F20.10)',('LA ',COORDSA(3*(J1-1)+1),COORDSA(3*(J1-1)+2),COORDSA(3*(J1-1)+3),J1=1,NATOMS)
! PRINT '(I6)',NATOMS
! PRINT '(A)','coordsb in minpermdist:'
! PRINT '(A,3F20.10)',('LB ',COORDSB(3*(J1-1)+1),COORDSB(3*(J1-1)+2),COORDSB(3*(J1-1)+3),J1=1,NATOMS)

DISTANCE=SQRT(DISTANCE) ! now changed to return distance, not distance^2 22/11/10 DJW

RETURN
END SUBROUTINE MINPERMDIST
