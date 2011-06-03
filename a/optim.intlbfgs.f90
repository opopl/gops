!   Copyright (C) 2003-2010 David J. Wales
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
SUBROUTINE INTLBFGS(QSTART,QFINISH,NMINFOUND,NTSFOUND,MINFOUND,TSFOUND)
USE PORFUNCS
USE KEY, ONLY : FREEZENODEST, FREEZETOL, DEBUG, MAXINTBFGS, NNREPULSIVE, &
     & INTCONSEP, INTRMSTOL, INTIMAGE, NREPMAX, NREPULSIVE, INTMUPDATE, INTDGUESS, &
     & INTCONSTRAINTT, NCONSTRAINT, CONI, CONJ, CONDISTREF, INTCONMAX, BSMIN, &
     & INTCONSTRAINREPCUT, REPCON, INTCONSTRAINTREP, INTREPSEP, NREPI, NREPJ, &
     & INTCONSTRAINTTOL, CONDISTREFLOCAL, INTCONFRAC, CONACTIVE, NITSTART, REPI, &
     & REPJ, NREPMAX, ATOMACTIVE, NCONSTRAINTON, CONION, CONJON, CONDISTREFLOCALON, CONDISTREFON, &
     & MAXCONUSE, NREPCUT, REPCUT, CHECKCONINT, INTCONSTEPS, INTRELSTEPS, MAXCONE, COLDFUSIONLIMIT, &
     & INTSTEPS1, DUMPINTXYZ, DUMPINTXYZFREQ, DUMPINTEOS, DUMPINTEOSFREQ, MUPDATE, BFGSSTEPS, INTTST, &
     & BFGSTST, NSTEPS, IMSEPMIN, IMSEPMAX, MAXINTIMAGE
USE COMMONS, ONLY: NATOMS, NOPT, ZSYM
USE MODEFOL

IMPLICIT NONE 

DOUBLE PRECISION, INTENT(IN) :: QSTART(NOPT), QFINISH(NOPT)  ! The two end points
INTEGER D, U
DOUBLE PRECISION DMAX, DS, DF, FRAC, DMIN
DOUBLE PRECISION, ALLOCATABLE :: REPTEMP(:)
INTEGER, ALLOCATABLE :: IREPTEMP(:)
INTEGER NDECREASE, NFAIL, NTOADD, NADDED, NMAXINT, NMININT, JMAX, JMIN, INTIMAGESAVE
LOGICAL KNOWE, KNOWG, KNOWH, ADDATOM, PTEST, MFLAG, PRINTOPTIMIZETS
COMMON /KNOWN/ KNOWE, KNOWG, KNOWH

DOUBLE PRECISION EDUMMY,EVALMIN,EVALMAX
INTEGER J2,POINT,BOUND,NPT,CP,I,ISTAT,J1,J3,J4,NIMAGEFREEZE,NACTIVE,NBEST,NEWATOM,NMINFOUND,NTSFOUND,NSIDE,ITDONE
INTEGER NCONFORNEWATOM, CONLIST(NATOMS), NCONATOM(NATOMS), TURNONORDER(NATOMS),NBACKTRACK
DOUBLE PRECISION, DIMENSION(3*NATOMS) :: LGDUMMY, VECS, XDIAG
INTEGER NDUMMY, CLOSEATOM, NLASTGOODE, NCONTOACTIVE(NATOMS), NSTEPSMAX
INTEGER NTRIES(NATOMS), N1, N2, N3, NITERDONE, EXITSTATUS
DOUBLE PRECISION :: YS,YY,SQ,YR,BETA,GNORM,DDOT,DUMMY,STPMIN,PREVGRAD,EMINUS,EPLUS, STARTTIME, TIME0, &
  &                 C1, C2, C3, VEC1(3), VEC2(3), VEC3(3), &
  &                 DUMMY2, RANDOM, INVDTOACTIVE(NATOMS), DINCREMENT, &
  &                 CONDIST(NATOMS), ESAVED, ESAVEC, ESAVE0, &
  &                 ETOTALTMP, RMSTMP, USEFRAC, & 
  &                 RAN1, STIME, FTIME, MAXCONDIST, &
  &                 MINCONDIST, ETOTAL, LASTGOODE, DPRAND, RMS, STEPTOT, &
  &                 VNEW(NOPT), ENERGY, RMS2, EREAL, MINCOORDS(2,NOPT), LINTCONSTRAINTTOL, LOCALCOORDS(3*NATOMS)

LOGICAL TSCONVERGED
DOUBLE PRECISION, DIMENSION(INTMUPDATE)     :: RHO1,ALPHA
DOUBLE PRECISION :: EOLD, DIFF, DIST, DTOTAL
LOGICAL SWITCHED, CHANGED
INTEGER NDIST1(NATOMS), NCYCLE, DMIN1, DMAX1, NUNCON1
DOUBLE PRECISION, POINTER :: X(:), G(:)
DOUBLE PRECISION, ALLOCATABLE :: EINT(:)
!
! These declarations have to match those in NEB/ntc.f90
!
TYPE MINFOUNDTYPE
   DOUBLE PRECISION,POINTER :: E
   DOUBLE PRECISION,POINTER :: COORD(:)
END TYPE MINFOUNDTYPE
INTEGER,PARAMETER :: NMINMAX = 3000 ! Maximal number of min to be checked in one intlbfgs run
TYPE (MINFOUNDTYPE) :: MINFOUND(NMINMAX)

INTEGER,PARAMETER :: NTSMAX = 3000 ! Maximal number of ts to be checked in one intlbfgs run
TYPE TSFOUNDTYPE
     DOUBLE PRECISION,POINTER :: E
     DOUBLE PRECISION,POINTER :: EVALMIN
     DOUBLE PRECISION,POINTER :: COORD(:)
     DOUBLE PRECISION,POINTER :: VECS(:)
END TYPE TSFOUNDTYPE

TYPE (TSFOUNDTYPE) :: TSFOUND(NTSMAX)
!
! efk: for freezenodes
!
DOUBLE PRECISION :: TESTG, TOTGNORM
INTEGER :: IM
INTEGER NDFORNEWATOM, BESTPRESERVEDN(NATOMS)
INTEGER NCFORNEWATOM, BESTCLOSESTN(NATOMS)
DOUBLE PRECISION BESTPRESERVEDD(NATOMS), BESTCLOSESTD(NATOMS)
!
! Dimensions involving INTIMAGE
!
DOUBLE PRECISION, ALLOCATABLE :: SAVEX(:), TRUEEE(:), GOODSAVE(:), &
  &              EEETMP(:), MYGTMP(:), XSAVED(:,:), XSAVEC(:,:), EEE(:), STEPIMAGE(:), &
  &              XSAVE0(:,:), GTMP(:), DIAG(:), STP(:), SEARCHSTEP(:,:), GDIF(:,:), GLAST(:), XSAVE(:)
DOUBLE PRECISION, ALLOCATABLE, TARGET :: XYZ(:), GGG(:), DPTMP(:), D2TMP(:,:)
LOGICAL, ALLOCATABLE :: CHECKG(:), IMGFREEZE(:), LOGTMP(:)

ALLOCATE(SAVEX(3*NATOMS*INTIMAGE),TRUEEE(INTIMAGE+2),GOODSAVE(3*NATOMS*(INTIMAGE+2)), &
  &      EEETMP(INTIMAGE+2), MYGTMP(3*NATOMS*INTIMAGE), XSAVED(3,INTIMAGE+2), &
  &      XSAVEC(3,INTIMAGE+2), XSAVE0(3,INTIMAGE+2), GTMP(3*NATOMS*INTIMAGE), &
  &      DIAG(3*NATOMS*INTIMAGE), STP(3*NATOMS*INTIMAGE), SEARCHSTEP(0:INTMUPDATE,NOPT*INTIMAGE), &
  &      GDIF(0:INTMUPDATE,NOPT*INTIMAGE),GLAST(NOPT*INTIMAGE), XSAVE(NOPT*INTIMAGE), &
  &      XYZ(NOPT*(INTIMAGE+2)), GGG(NOPT*(INTIMAGE+2)), CHECKG(NOPT*INTIMAGE), IMGFREEZE(INTIMAGE), &
  &      EEE(INTIMAGE+2), STEPIMAGE(INTIMAGE))

SWITCHED=.FALSE.
INTIMAGESAVE=INTIMAGE
NBACKTRACK=1
CALL MYCPU_TIME(STIME,.FALSE.)
PRINT '(A,I6)',' intlbfgs> Maximum number of steps for constraint potential phase is ',INTSTEPS1
PREVGRAD=1.0D100
ADDATOM=.FALSE.
NFAIL=0
IF (FREEZENODEST) IMGFREEZE(1:INTIMAGE)=.FALSE.
D=NOPT*INTIMAGE
U=INTMUPDATE
NITERDONE=1
LINTCONSTRAINTTOL=INTCONSTRAINTTOL

IF ( D<=0 ) THEN
   PRINT *, 'd is not positive, d=',d
   CALL TSUMMARY
   STOP
ENDIF
IF ( U<=0 ) THEN
   PRINT *, 'u is not positive, u=',u
   CALL TSUMMARY
   STOP
ENDIF
IF (INTSTEPS1 < 0) THEN
   PRINT '(1x,a)', 'Maximal number of iterations is less than zero! Stop.'
   CALL TSUMMARY
   STOP
ENDIF
!
! XYZ, GGG, EEE include the end point images
! X, G do not.
!
IF (.NOT.ALLOCATED(CONI)) THEN 
   ALLOCATE(CONI(INTCONMAX),CONJ(INTCONMAX),CONDISTREF(INTCONMAX))
   ALLOCATE(REPI(NREPMAX),REPJ(NREPMAX),NREPI(NREPMAX),NREPJ(NREPMAX),REPCUT(NREPMAX),NREPCUT(NREPMAX))
ENDIF
X=>XYZ(NOPT+1:NOPT*(INTIMAGE+1))
G=>GGG(NOPT+1:NOPT*(INTIMAGE+1))
!
! Initialise XYZ
!
XYZ(1:NOPT)=QSTART(1:NOPT)
XYZ(NOPT*(INTIMAGE+1)+1:NOPT*(INTIMAGE+2))=QFINISH(1:NOPT)
DO J1=1,INTIMAGE+2
   XYZ((J1-1)*NOPT+1:J1*NOPT)=((INTIMAGE+2-J1)*QSTART(1:NOPT)+(J1-1)*QFINISH(1:NOPT))/(INTIMAGE+1)
ENDDO
!
! Calculate initial constraints.
!
NLASTGOODE=0
LASTGOODE=1.0D100
GOODSAVE(1:D)=X(1:D)
IF (.NOT.ALLOCATED(ATOMACTIVE)) ALLOCATE(ATOMACTIVE(NATOMS))
51   NCONSTRAINT=0 
MAXCONDIST=-1.0D0
SAVEX(1:D)=X(1:D)
MINCONDIST=1.0D100
DO J2=1,NATOMS
   DO J3=J2+1,NATOMS
      IF (J3-J2.GT.INTCONSEP) CYCLE ! forbid constraints corresponding to atoms distant in sequence
      DS=SQRT((XYZ(3*(J2-1)+1)-XYZ(3*(J3-1)+1))**2 &
  &          +(XYZ(3*(J2-1)+2)-XYZ(3*(J3-1)+2))**2 &
  &          +(XYZ(3*(J2-1)+3)-XYZ(3*(J3-1)+3))**2) 
      IF (DS.GT.5.0D0) CYCLE ! don't allow constraints if either endpoint separation is too large DJW
      DF=SQRT((XYZ((INTIMAGE+1)*3*NATOMS+3*(J2-1)+1)-XYZ((INTIMAGE+1)*3*NATOMS+3*(J3-1)+1))**2 &
  &          +(XYZ((INTIMAGE+1)*3*NATOMS+3*(J2-1)+2)-XYZ((INTIMAGE+1)*3*NATOMS+3*(J3-1)+2))**2 &
  &          +(XYZ((INTIMAGE+1)*3*NATOMS+3*(J2-1)+3)-XYZ((INTIMAGE+1)*3*NATOMS+3*(J3-1)+3))**2) 
      IF (DF.GT.5.0D0) CYCLE ! don't allow constraints if either endpoint separation is too large DJW
!     IF (2.0D0*ABS(DS-DF)/(DS+DF).LT.LINTCONSTRAINTTOL) THEN
      IF (ABS(DS-DF).LT.LINTCONSTRAINTTOL) THEN
!
!  Add constraint for this distance to the list.
!
         NCONSTRAINT=NCONSTRAINT+1
!        PRINT '(A,2I6,A,I6)','intlbfgs> Adding constraint for atoms ',J2,J3,'  total=',NCONSTRAINT
         IF (NCONSTRAINT.GT.INTCONMAX) THEN
            ALLOCATE(IREPTEMP(INTCONMAX))
               
            IREPTEMP(1:INTCONMAX)=CONI(1:INTCONMAX)
            DEALLOCATE(CONI)
            ALLOCATE(CONI(2*INTCONMAX))
            CONI(1:INTCONMAX)=IREPTEMP(1:INTCONMAX)
               
            IREPTEMP(1:INTCONMAX)=CONJ(1:INTCONMAX)
            DEALLOCATE(CONJ)
            ALLOCATE(CONJ(2*INTCONMAX))
            CONJ(1:INTCONMAX)=IREPTEMP(1:INTCONMAX)
               
            DEALLOCATE(IREPTEMP)
            ALLOCATE(REPTEMP(1:INTCONMAX))
               
            REPTEMP(1:INTCONMAX)=CONDISTREF(1:INTCONMAX)
            DEALLOCATE(CONDISTREF)
            ALLOCATE(CONDISTREF(2*INTCONMAX))
            CONDISTREF(1:INTCONMAX)=REPTEMP(1:INTCONMAX)

            INTCONMAX=2*INTCONMAX
            DEALLOCATE(REPTEMP)
         ENDIF
         CONI(NCONSTRAINT)=J2
         CONJ(NCONSTRAINT)=J3
         CONDISTREF(NCONSTRAINT)=(DF+DS)/2.0D0
         IF (CONDISTREF(NCONSTRAINT).GT.MAXCONDIST) MAXCONDIST=CONDISTREF(NCONSTRAINT)
         IF (CONDISTREF(NCONSTRAINT).LT.MINCONDIST) MINCONDIST=CONDISTREF(NCONSTRAINT)
         IF (DEBUG) PRINT '(A,2I6,A,2F12.2,A,F12.4,A,I8)',' intlbfgs> constrain distance for atoms ',CONI(NCONSTRAINT), &
  &              CONJ(NCONSTRAINT),' values are ',DS,DF,' fraction=',2*ABS(DS-DF)/(DS+DF), &
  &             ' # constraints=',NCONSTRAINT
      ENDIF
   ENDDO
ENDDO
!
! Check that we have a percolating constraint network. If not, increase the tolerance and try again!
! Calculate minimum number of steps of each atom from number 1.
!
NDIST1(1:NATOMS)=1000000
NDIST1(1)=0
NCYCLE=0
5    CHANGED=.FALSE.
NCYCLE=NCYCLE+1
DMIN1=100000
DMAX1=0
NUNCON1=0
DO J1=1,NATOMS
   IF (NDIST1(J1).EQ.0) CYCLE ! minimum 1
     DO J2=1,NCONSTRAINT
        IF (CONI(J2).EQ.J1) THEN
           IF (NDIST1(CONJ(J2))+1.LT.NDIST1(J1)) THEN
              CHANGED=.TRUE.
              NDIST1(J1)=NDIST1(CONJ(J2))+1
           ENDIF
        ELSE IF (CONJ(J2).EQ.J1) THEN
           IF (NDIST1(CONI(J2))+1.LT.NDIST1(J1)) THEN
              CHANGED=.TRUE.
              NDIST1(J1)=NDIST1(CONI(J2))+1
           ENDIF
        ENDIF
     ENDDO
     IF ((NDIST1(J1).GT.DMAX1).AND.(NDIST1(J1).NE.1000000)) DMAX1=NDIST1(J1)
     IF (NDIST1(J1).LT.DMIN1) DMIN1=NDIST1(J1)
     IF (NDIST1(J1).EQ.1000000) NUNCON1=NUNCON1+1
ENDDO
IF (CHANGED) GOTO 5
PRINT '(3(A,I8))',' intlbfgs> steps to atom 1 converged in ',NCYCLE-1, &
  &               ' cycles; maximum=',DMAX1,' disconnected=',NUNCON1
IF (NUNCON1.GT.0) THEN
   LINTCONSTRAINTTOL=LINTCONSTRAINTTOL*1.1D0
   PRINT '(A,F15.5)',' intlbfgs> increasing the local constraint tolerance parameter to ',LINTCONSTRAINTTOL
   GOTO 51
ENDIF
NCONATOM(1:NATOMS)=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PRINT '(A,I6,2(A,F15.5))',' intlbfgs> total distance constraints=',NCONSTRAINT, &
  &                    ' shortest=',MINCONDIST,' longest=',MAXCONDIST
REPCON=-INTCONSTRAINTREP/INTCONSTRAINREPCUT**6 ! also needed for congrad.f90 potential
IF (ALLOCATED(CONDISTREFLOCAL)) THEN
   DEALLOCATE(CONDISTREFLOCAL)
ENDIF
ALLOCATE(CONDISTREFLOCAL(NCONSTRAINT))
IF (ALLOCATED(CONDISTREFLOCALON)) DEALLOCATE(CONDISTREFLOCALON)
IF (ALLOCATED(CONDISTREFON)) DEALLOCATE(CONDISTREFON)
IF (ALLOCATED(CONION)) DEALLOCATE(CONION)
IF (ALLOCATED(CONJON)) DEALLOCATE(CONJON)
ALLOCATE(CONDISTREFLOCALON(NCONSTRAINT),CONDISTREFON(NCONSTRAINT),CONION(NCONSTRAINT),CONJON(NCONSTRAINT))
CONDISTREFLOCAL(1:NCONSTRAINT)=CONDISTREF(1:NCONSTRAINT)
DUMMY=1.0D100
DO J1=1,NCONSTRAINT
   DF=SQRT((XYZ(3*(CONI(J1)-1)+1)-XYZ((INTIMAGE+1)*3*NATOMS+3*(CONI(J1)-1)+1))**2 &
  &       +(XYZ(3*(CONI(J1)-1)+2)-XYZ((INTIMAGE+1)*3*NATOMS+3*(CONI(J1)-1)+2))**2 &
  &       +(XYZ(3*(CONI(J1)-1)+3)-XYZ((INTIMAGE+1)*3*NATOMS+3*(CONI(J1)-1)+3))**2)&
  &  +SQRT((XYZ(3*(CONJ(J1)-1)+1)-XYZ((INTIMAGE+1)*3*NATOMS+3*(CONJ(J1)-1)+1))**2 &
  &       +(XYZ(3*(CONJ(J1)-1)+2)-XYZ((INTIMAGE+1)*3*NATOMS+3*(CONJ(J1)-1)+2))**2 &
  &       +(XYZ(3*(CONJ(J1)-1)+3)-XYZ((INTIMAGE+1)*3*NATOMS+3*(CONJ(J1)-1)+3))**2)
   IF (DF.LT.DUMMY) THEN
      NBEST=J1
      DUMMY=DF
   ENDIF
ENDDO
IF (DEBUG) PRINT '(A,I6,A,2I6,A,F15.5)',' intlbfgs> Smallest overall motion for constraint ',NBEST,' atoms ', &
  &                           CONI(NBEST),CONJ(NBEST),' distance=',DUMMY
NACTIVE=2
TURNONORDER(1:NATOMS)=0
NTRIES(1:NATOMS)=1
IF (ALLOCATED(CONACTIVE)) DEALLOCATE(CONACTIVE)
IF (ALLOCATED(NITSTART)) DEALLOCATE(NITSTART)
ALLOCATE(CONACTIVE(NCONSTRAINT),NITSTART(NCONSTRAINT))
CONACTIVE(1:NCONSTRAINT)=.FALSE.
CONACTIVE(NBEST)=.TRUE.
NITSTART(NBEST)=1
ATOMACTIVE(1:NATOMS)=.FALSE.
ATOMACTIVE(CONI(NBEST))=.TRUE.
ATOMACTIVE(CONJ(NBEST))=.TRUE.
TURNONORDER(1)=CONI(NBEST)
TURNONORDER(2)=CONJ(NBEST)
NTRIES(CONI(NBEST))=1
NTRIES(CONJ(NBEST))=1
NREPULSIVE=0
NCONSTRAINTON=1
CONDISTREFLOCALON(1)=CONDISTREFLOCAL(NBEST)
CONDISTREFON(1)=CONDISTREF(NBEST)
CONION(1)=CONI(NBEST)
CONJON(1)=CONJ(NBEST)
!
! Don;t want to redistribute images before even taking a step, so don;t call CHECKSEP.
! Must call CHECKREP to initialise NNREULSIVE, NREPI, NREPJ, etc. SEGV otherwise on second cycle!
!

CALL CHECKREP(INTIMAGE,XYZ,NOPT)
IF (CHECKCONINT) THEN
   CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
ELSE
   CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
ENDIF
EOLD=ETOTAL
GLAST(1:D)=G(1:D)
XSAVE(1:D)=X(1:D)

IF (ETOTAL/INTIMAGE.LT.COLDFUSIONLIMIT) THEN
   WRITE(*,'(A,2G20.10)') ' intlbfgs> Cold fusion diagnosed - step discarded, energy, limit=', &
  &                       ETOTAL/INTIMAGE,COLDFUSIONLIMIT
   DEALLOCATE(CONI,CONJ,CONDISTREF,REPI,REPJ,NREPI,NREPJ,REPCUT,NREPCUT)
   DEALLOCATE(SAVEX,TRUEEE,GOODSAVE, EEETMP, MYGTMP, XSAVED, XSAVEC, XSAVE0, GTMP, &
  &      DIAG, STP, SEARCHSTEP, GDIF,GLAST, XSAVE, XYZ, GGG, CHECKG, IMGFREEZE, EEE, STEPIMAGE)
   INTIMAGE=INTIMAGESAVE
   NTSFOUND=0
   NMINFOUND=0
   RETURN
ENDIF

IF (DEBUG) WRITE(*,'(A6,A20,A20,A9,A9)') 'Iter','Energy per image','RMS Force','Step'
NSTEPSMAX=INTSTEPS1


DO ! Main do loop with counter NITERDONE, initially set to one
!
!  Add next atom to active set if ADDATOM is true. 
!  Constraints to atoms already in the active set are turned on
!  and short-range repulsions to active atoms that are not distance constrained are turned on.
!  *** OLD Find nearest atom to active set attached by a constraint
!  *** NEW Find atom with most constraints to active set
!  Turn on constraint terms for this atom with all previous members of the active set
!  Add repulsions to non-constrained atoms in this set
!  NTOADD is the number of atoms to add to the active set in each pass. 1 seems best!
!
   NTOADD=1
!  NTOADD=NATOMS-2  !!!! DJW
   NADDED=0
   IF (ADDATOM.AND.(NACTIVE.LT.NATOMS)) THEN
542   CONTINUE
!     DUMMY=1.0D100
      NBEST=0
      NCONTOACTIVE(1:NATOMS)=0
      INVDTOACTIVE(1:NATOMS)=0.0D0
      SAVEX(D)=X(D)
      DO J2=1,NCONSTRAINT
         IF (CONACTIVE(J2)) CYCLE   ! count new, inactive constraints
         IF (ATOMACTIVE(CONI(J2))) THEN
            IF (.NOT.ATOMACTIVE(CONJ(J2))) THEN
               NCONTOACTIVE(CONJ(J2))=NCONTOACTIVE(CONJ(J2))+1
               INVDTOACTIVE(CONJ(J2))=INVDTOACTIVE(CONJ(J2))+1.0D0/CONDISTREF(J2)
            ENDIF
         ENDIF
         IF (ATOMACTIVE(CONJ(J2))) THEN
            IF (.NOT.ATOMACTIVE(CONI(J2))) THEN
               NCONTOACTIVE(CONI(J2))=NCONTOACTIVE(CONI(J2))+1
               INVDTOACTIVE(CONI(J2))=INVDTOACTIVE(CONI(J2))+1.0D0/CONDISTREF(J2)
            ENDIF
         ENDIF
         IF (NCONTOACTIVE(CONI(J2)).GT.NBEST) THEN
            NBEST=NCONTOACTIVE(CONI(J2))
         ENDIF
         IF (NCONTOACTIVE(CONJ(J2)).GT.NBEST) THEN
            NBEST=NCONTOACTIVE(CONJ(J2))
         ENDIF
!        PRINT '(A,7I6)','J2,NCONTOACTIVEI,NCONTOACTOVEJ,CONI,CONJ,NEWATOM,NBEST=', &
! &                             J2,NCONTOACTIVE(CONI(J2)),NCONTOACTIVE(CONJ(J2)),CONI(J2),CONJ(J2),NEWATOM,NBEST

      ENDDO
!
!  Choose NEWATOM stochastically. Bias towards atoms with the maximum constraints.
!  Use a normalised probability and generate a random number between 0 and 1.
!
      DUMMY2=0.0D0
      DO J2=1,NATOMS
         IF (NCONTOACTIVE(J2).EQ.0) CYCLE
         IF (ATOMACTIVE(J2)) CYCLE
!        DUMMY2=DUMMY2+((1.0D0*NCONTOACTIVE(J2))/(1.0D0*CONDISTREF(J2)*NTRIES(J2)))**4 
         DUMMY2=DUMMY2+((1.0D0*INVDTOACTIVE(J2))/(1.0D0*NTRIES(J2)))**4 
!        PRINT '(A,I6,A,G20.10)',' intlbfgs> Unnormalised probability for choosing atom ',J2,' is ', &
! &                ((1.0D0*INVDTOACTIVE(J2))/(1.0D0*NTRIES(J2)))**4
      ENDDO

      RANDOM=DUMMY2*DPRAND()
      DUMMY2=0.0D0
      choosenew: DO J2=1,NATOMS
         IF (NCONTOACTIVE(J2).EQ.0) CYCLE
         IF (ATOMACTIVE(J2)) CYCLE
!        DUMMY2=DUMMY2+((1.0D0*NCONTOACTIVE(J2))/(1.0D0*CONDISTREF(J2)*NTRIES(J2)))**4 
         DUMMY2=DUMMY2+((1.0D0*INVDTOACTIVE(J2))/(1.0D0*NTRIES(J2)))**4 
         IF (DUMMY2.GE.RANDOM) THEN
            NEWATOM=J2
            IF (DEBUG) PRINT '(3(A,I6))',' intlbfgs> Choosing new active atom ',NEWATOM,' new constraints=', &
  &                                       NCONTOACTIVE(J2),' maximum=',NBEST
            EXIT choosenew
         ENDIF
      ENDDO choosenew
          
      IF (NEWATOM*NBEST.EQ.0) THEN ! sanity check
         PRINT '(A,I6,A,2I6)',' intlbfgs> ERROR *** new active atom not set'
         STOP
      ELSE
!
!  We need a sorted list of up to 3 active atoms, sorted according to how well the
!  end point distance is preserved, even if they don't satisfy the constraint 
!  condition. We want three atoms to use for a local axis system in the interpolation.
!
!  Try sorting on the shortest average distances in the endpoint structures instead, to avoid
!  problems with distant atoms acidentally having a well-preserved distance.
!
         NDFORNEWATOM=0
         BESTPRESERVEDD(1:NATOMS)=1.0D100
         DO J1=1,NATOMS
            IF (.NOT.ATOMACTIVE(J1)) CYCLE
            DS=SQRT((XYZ(3*(NEWATOM-1)+1)-XYZ(3*(J1-1)+1))**2 &
  &                +(XYZ(3*(NEWATOM-1)+2)-XYZ(3*(J1-1)+2))**2 &
  &                +(XYZ(3*(NEWATOM-1)+3)-XYZ(3*(J1-1)+3))**2) 
            DF=SQRT((XYZ((INTIMAGE+1)*3*NATOMS+3*(NEWATOM-1)+1)-XYZ((INTIMAGE+1)*3*NATOMS+3*(J1-1)+1))**2 &
  &                +(XYZ((INTIMAGE+1)*3*NATOMS+3*(NEWATOM-1)+2)-XYZ((INTIMAGE+1)*3*NATOMS+3*(J1-1)+2))**2 &
  &                +(XYZ((INTIMAGE+1)*3*NATOMS+3*(NEWATOM-1)+3)-XYZ((INTIMAGE+1)*3*NATOMS+3*(J1-1)+3))**2) 
            DUMMY=ABS(DS-DF)
            NDFORNEWATOM=NDFORNEWATOM+1
            DO J2=1,NDFORNEWATOM 
               IF (DUMMY.LT.BESTPRESERVEDD(J2)) THEN
!                 PRINT '(A,I6,G12.4,I6,G12.4)','J1,DUMMY < J2,BESTPRESERVEDD: ',J1,DUMMY,J2,BESTPRESERVEDD(J2)
                  DO J3=NDFORNEWATOM,J2+1,-1 
!                    PRINT '(A,I6,A,I6,A,G12.4)',' moving diff and list from ',J3-1,' to ',J3, &
!&                                               ' DIFF=',BESTPRESERVEDD(J3-1)
                     BESTPRESERVEDD(J3)=BESTPRESERVEDD(J3-1)
                     BESTPRESERVEDN(J3)=BESTPRESERVEDN(J3-1)
                  ENDDO
                  BESTPRESERVEDD(J2)=DUMMY
!                 PRINT '(A,I6,A,G12.4)',' setting BESTPRESERVEDD element ',J2,' to ',DUMMY
                  BESTPRESERVEDN(J2)=J1
!                 PRINT '(A,I6,A,G12.4)',' setting BESTPRESERVEDN element ',J2,' to ',J1
                  GOTO 653
               ENDIF
            ENDDO
653         CONTINUE
         ENDDO
         IF (DEBUG) THEN
            PRINT '(A,I6,A,I6,A)',' intlbfgs> New active atom ',NEWATOM,' best preserved distances:'
            PRINT '(20I6)',BESTPRESERVEDN(1:MIN(10,NDFORNEWATOM))
            PRINT '(A,I6,A,I6,A)',' intlbfgs> sorted differences:'
            PRINT '(10G12.4)',BESTPRESERVEDD(1:MIN(10,NDFORNEWATOM))
         ENDIF
         IF (FREEZENODEST) IMGFREEZE(1:INTIMAGE)=.FALSE.

         NCFORNEWATOM=0
         BESTCLOSESTD(1:NATOMS)=1.0D100
         DO J1=1,NATOMS
            IF (.NOT.ATOMACTIVE(J1)) CYCLE
            DS=SQRT((XYZ(3*(NEWATOM-1)+1)-XYZ(3*(J1-1)+1))**2 &
  &                +(XYZ(3*(NEWATOM-1)+2)-XYZ(3*(J1-1)+2))**2 &
  &                +(XYZ(3*(NEWATOM-1)+3)-XYZ(3*(J1-1)+3))**2) 
            DF=SQRT((XYZ((INTIMAGE+1)*3*NATOMS+3*(NEWATOM-1)+1)-XYZ((INTIMAGE+1)*3*NATOMS+3*(J1-1)+1))**2 &
  &                +(XYZ((INTIMAGE+1)*3*NATOMS+3*(NEWATOM-1)+2)-XYZ((INTIMAGE+1)*3*NATOMS+3*(J1-1)+2))**2 &
  &                +(XYZ((INTIMAGE+1)*3*NATOMS+3*(NEWATOM-1)+3)-XYZ((INTIMAGE+1)*3*NATOMS+3*(J1-1)+3))**2) 
            DUMMY=(DS+DF)/2.0D0
            NCFORNEWATOM=NCFORNEWATOM+1
            DO J2=1,NCFORNEWATOM
               IF (DUMMY.LT.BESTCLOSESTD(J2)) THEN
!                 PRINT '(A,I6,G12.4,I6,G12.4)','J1,DUMMY < J2,BESTCLOSESTD: ',J1,DUMMY,J2,BESTCLOSESTD(J2)
                  DO J3=NCFORNEWATOM,J2+1,-1
!                    PRINT '(A,I6,A,I6,A,G12.4)',' moving diff and list from ',J3-1,' to ',J3, &
!&                                               ' DIFF=',BESTCLOSESTD(J3-1)
                     BESTCLOSESTD(J3)=BESTCLOSESTD(J3-1)
                     BESTCLOSESTN(J3)=BESTCLOSESTN(J3-1)
                  ENDDO
                  BESTCLOSESTD(J2)=DUMMY
!                 PRINT '(A,I6,A,G12.4)',' setting BESTCLOSESTD element ',J2,' to ',DUMMY
                  BESTCLOSESTN(J2)=J1
!                 PRINT '(A,I6,A,G12.4)',' setting BESTCLOSESTN element ',J2,' to ',J1
                  GOTO 659
               ENDIF
            ENDDO
659         CONTINUE
         ENDDO
         IF (DEBUG) THEN
            PRINT '(A,I6,A,I6,A)',' intlbfgs> New active atom ',NEWATOM,' shortest average distances in endpoints:'
            PRINT '(20I6)',BESTCLOSESTN(1:MIN(10,NCFORNEWATOM))
            PRINT '(A,I6,A,I6,A)',' intlbfgs> sorted differences:'
            PRINT '(10G12.4)',BESTCLOSESTN(1:MIN(10,NCFORNEWATOM))
         ENDIF
!
!  Maintain a sorted list of active atoms that are constrained to the new atom, sorted
!  according to their distance.
!
         NCONFORNEWATOM=0
         CONDIST(1:NATOMS)=1.0D100
         IF (DEBUG) PRINT '(3(A,I6))',' intlbfgs> New active atom is number ',NEWATOM,' total=',NACTIVE+1, &
 &                        ' steps=',NITERDONE
         DO J1=1,NCONSTRAINT
            IF (CONACTIVE(J1)) CYCLE
            IF ((CONI(J1).EQ.NEWATOM).AND.(ATOMACTIVE(CONJ(J1))).OR.(CONJ(J1).EQ.NEWATOM).AND.(ATOMACTIVE(CONI(J1)))) THEN  
                 NCONFORNEWATOM=NCONFORNEWATOM+1
!                CONACTIVE(J1)=.TRUE.
!                NITSTART(J1)=NITERDONE
!                NCONSTRAINTON=NCONSTRAINTON+1
! !
! ! The ...ON variables are not actually used in congrad.f90.
! !
!                CONDISTREFLOCALON(NCONSTRAINTON)=CONDISTREFLOCAL(J1)
!                CONDISTREFON(NCONSTRAINTON)=CONDISTREF(J1)
!                CONION(NCONSTRAINTON)=CONI(J1)
!                CONJON(NCONSTRAINTON)=CONJ(J1)
! 
!                IF (DEBUG) PRINT '(A,I6,A,2I6)',' intlbfgs> Turning on constraint ',J1,' for atoms ',CONI(J1),CONJ(J1)
               IF (NCONFORNEWATOM.EQ.1) THEN
                  CONDIST(1)=CONDISTREF(J1)
                  IF (CONI(J1).EQ.NEWATOM) CONLIST(1)=CONJ(J1)
                  IF (CONJ(J1).EQ.NEWATOM) CONLIST(1)=CONI(J1)
               ENDIF
               DO J2=1,NCONFORNEWATOM-1
                  IF (CONDISTREF(J1).LT.CONDIST(J2)) THEN
!                    PRINT '(A,I6,G12.4,I6,G12.4)','J1,CONDISTREF < J2,CONDIST: ',J1,CONDISTREF(J1),J2,CONDIST(J2)
                     DO J3=NCONFORNEWATOM,J2+1,-1
!                       PRINT '(A,I6,A,I6,A,G12.4)',' moving dist and list from ',J3-1,' to ',J3,' CONDIST=',CONDIST(J3-1)
                        CONDIST(J3)=CONDIST(J3-1)
                        CONLIST(J3)=CONLIST(J3-1)
                     ENDDO
                     CONDIST(J2)=CONDISTREF(J1)
!                    PRINT '(A,I6,A,G12.4)',' setting condist element ',J2,' to ',CONDISTREF(J1)
                     IF (CONI(J1).EQ.NEWATOM) CONLIST(J2)=CONJ(J1)
                     IF (CONJ(J1).EQ.NEWATOM) CONLIST(J2)=CONI(J1)
!                    PRINT '(A,I6,A,G12.4)',' setting conlist element ',J2,' to ',CONLIST(J2)
                     GOTO 654
                  ENDIF
               ENDDO 
               CONDIST(NCONFORNEWATOM)=CONDISTREF(J1)
!              PRINT '(A,I6,A,G12.4)',' setting condist element ',NCONFORNEWATOM,' to ',CONDISTREF(J1)
               IF (CONI(J1).EQ.NEWATOM) CONLIST(NCONFORNEWATOM)=CONJ(J1)
               IF (CONJ(J1).EQ.NEWATOM) CONLIST(NCONFORNEWATOM)=CONI(J1)
!              PRINT '(A,I6,A,G12.4)',' setting conlist element ',NCONFORNEWATOM,' to ',CONLIST(NCONFORNEWATOM)
654          CONTINUE
            ENDIF
         ENDDO 
         IF (DEBUG) THEN
            PRINT '(A,I6,A,I6,A)',' intlbfgs> New active atom ',NEWATOM,' is constrained to ',NCONFORNEWATOM,' other active atoms:'
            PRINT '(20I6)',CONLIST(1:NCONFORNEWATOM)
            PRINT '(A,I6,A,I6,A)',' intlbfgs> sorted distances:'
            PRINT '(10G12.4)',CONDIST(1:NCONFORNEWATOM)
         ENDIF
         DO J1=1,MIN(MAXCONUSE,NCONFORNEWATOM)
            DO J2=1,NCONSTRAINT
               IF ((CONI(J2).EQ.NEWATOM).AND.(CONJ(J2).EQ.CONLIST(J1))) THEN
                     CONACTIVE(J2)=.TRUE.
                     NCONATOM(CONI(J2))=NCONATOM(CONI(J2))+1
                     NCONATOM(CONJ(J2))=NCONATOM(CONJ(J2))+1
                     IF (DEBUG) PRINT '(A,I6,A,2I6)',' intlbfgs> Turning on constraint ',J2,' for atoms ',CONI(J2),CONJ(J2)
               ELSE IF ((CONJ(J2).EQ.NEWATOM).AND.(CONI(J2).EQ.CONLIST(J1))) THEN
                     CONACTIVE(J2)=.TRUE.
                     NCONATOM(CONI(J2))=NCONATOM(CONI(J2))+1
                     NCONATOM(CONJ(J2))=NCONATOM(CONJ(J2))+1
                     IF (DEBUG) PRINT '(A,I6,A,2I6)',' intlbfgs> Turning on constraint ',J2,' for atoms ',CONI(J2),CONJ(J2)
               ENDIF
            ENDDO
         ENDDO
         DO J1=1,NATOMS
            IF (.NOT.ATOMACTIVE(J1)) CYCLE ! identify active atoms
            IF (ABS(J1-NEWATOM).LE.INTREPSEP) CYCLE ! no repulsion for atoms too close in sequence
            DO J2=1,NCONSTRAINT
!
!  With MAXCONUSE set to a finite value there could be constraints for the new atom that are
!  not active. We don't want these to be changed to repulsion, surely?!
!  Or perhaps we do need to do something with them?
!
               IF (.NOT.CONACTIVE(J2)) CYCLE ! identify active constraints 
               IF (((CONI(J2).EQ.J1).AND.(CONJ(J2).EQ.NEWATOM)).OR.((CONJ(J2).EQ.J1).AND.(CONI(J2).EQ.NEWATOM))) GOTO 543
            ENDDO
            DMIN=1.0D100
            DMAX=-1.0D0
            DO J2=1,INTIMAGE+2,INTIMAGE+1 ! only consider the end-point distances
               DF=SQRT((XYZ((J2-1)*3*NATOMS+3*(NEWATOM-1)+1)-XYZ((J2-1)*3*NATOMS+3*(J1-1)+1))**2+ &
  &                    (XYZ((J2-1)*3*NATOMS+3*(NEWATOM-1)+2)-XYZ((J2-1)*3*NATOMS+3*(J1-1)+2))**2+ &
  &                    (XYZ((J2-1)*3*NATOMS+3*(NEWATOM-1)+3)-XYZ((J2-1)*3*NATOMS+3*(J1-1)+3))**2)
               IF (DF.GT.DMAX) DMAX=DF
               IF (DF.LT.DMIN) DMIN=DF
            ENDDO
!
! Use the minimum of the end point distances and INTCONSTRAINREPCUT for each contact.
!
            DMIN=MIN(DMIN-1.0D-3,INTCONSTRAINREPCUT)
            NREPULSIVE=NREPULSIVE+1
            IF (NREPULSIVE.GT.NREPMAX) THEN
               ALLOCATE(IREPTEMP(NREPMAX),REPTEMP(NREPMAX))

               IREPTEMP(1:NREPMAX)=REPI(1:NREPMAX)
               DEALLOCATE(REPI)
               ALLOCATE(REPI(2*NREPMAX))
               REPI(1:NREPMAX)=IREPTEMP(1:NREPMAX)
               IREPTEMP(1:NREPMAX)=REPJ(1:NREPMAX)
               DEALLOCATE(REPJ)
               ALLOCATE(REPJ(2*NREPMAX))
               REPJ(1:NREPMAX)=IREPTEMP(1:NREPMAX)
               REPTEMP(1:NREPMAX)=REPCUT(1:NREPMAX)
               DEALLOCATE(REPCUT)
               ALLOCATE(REPCUT(2*NREPMAX))
               REPCUT(1:NREPMAX)=REPTEMP(1:NREPMAX)

               IREPTEMP(1:NREPMAX)=NREPI(1:NREPMAX)
               DEALLOCATE(NREPI)
               ALLOCATE(NREPI(2*NREPMAX))
               NREPI(1:NREPMAX)=IREPTEMP(1:NREPMAX)
               IREPTEMP(1:NREPMAX)=NREPJ(1:NREPMAX)
               DEALLOCATE(NREPJ)
               ALLOCATE(NREPJ(2*NREPMAX))
               NREPJ(1:NREPMAX)=IREPTEMP(1:NREPMAX)
               REPTEMP(1:NREPMAX)=NREPCUT(1:NREPMAX)
               DEALLOCATE(NREPCUT)
               ALLOCATE(NREPCUT(2*NREPMAX))
               NREPCUT(1:NREPMAX)=REPTEMP(1:NREPMAX)

               DEALLOCATE(IREPTEMP,REPTEMP)
               NREPMAX=2*NREPMAX
            ENDIF
            REPI(NREPULSIVE)=J1
            REPJ(NREPULSIVE)=NEWATOM
            REPCUT(NREPULSIVE)=DMIN
            IF (DEBUG) PRINT '(A,I6,A,I6,A,F15.5)',' intlbfgs> Adding repulsion for new atom ',NEWATOM,' with atom ',J1, &
  &                                                   ' cutoff=',DMIN
543         CONTINUE
         ENDDO
         ATOMACTIVE(NEWATOM)=.TRUE.
         NACTIVE=NACTIVE+1

         NDUMMY=0
         DO J1=1,NATOMS
            IF (ATOMACTIVE(J1)) NDUMMY=NDUMMY+1
         ENDDO
         IF (NDUMMY.NE.NACTIVE) THEN
            PRINT '(A,I6)',' intlbfgs> ERROR *** inconsistency in number of active atoms. Should be ',NACTIVE
            DO J1=1,NATOMS
               IF (ATOMACTIVE(J1)) PRINT '(A,I6)',' active atom ',J1
            ENDDO
            STOP
         ENDIF

         TURNONORDER(NACTIVE)=NEWATOM
!
! Initial guess for new active atom position. This is crucial for success in INTCONSTRAINT schemes!
!
         ESAVED=1.0D100
         ESAVE0=1.0D100
         ESAVEC=1.0D100
         IF (NCONFORNEWATOM.GE.3) THEN
!
! Move the new atom consistently in the local environment of its three nearest actively constrained atoms.
! Make a local orthogonal coordinate system and use constant components in this basis.
!
            IF (DEBUG) PRINT '(A)',' intlbfgs> initial guess from closest three constrained active atoms'
            VEC1(1:3)=XYZ(3*(CONLIST(2)-1)+1:3*(CONLIST(2)-1)+3)-XYZ(3*(CONLIST(1)-1)+1:3*(CONLIST(1)-1)+3)
            DUMMY=SQRT(VEC1(1)**2+VEC1(2)**2+VEC1(3)**2)
            IF (DUMMY.NE.0.0D0) VEC1(1:3)=VEC1(1:3)/DUMMY
            VEC2(1:3)=XYZ(3*(CONLIST(3)-1)+1:3*(CONLIST(3)-1)+3)-XYZ(3*(CONLIST(1)-1)+1:3*(CONLIST(1)-1)+3)
            DUMMY=VEC1(1)*VEC2(1)+VEC1(2)*VEC2(2)+VEC1(3)*VEC2(3)
            VEC2(1:3)=VEC2(1:3)-DUMMY*VEC1(1:3)
            DUMMY=SQRT(VEC2(1)**2+VEC2(2)**2+VEC2(3)**2)
            IF (DUMMY.NE.0.0D0) VEC2(1:3)=VEC2(1:3)/DUMMY
            VEC3(1)= VEC1(2)*VEC2(3)-VEC1(3)*VEC2(2)
            VEC3(2)=-VEC1(1)*VEC2(3)+VEC1(3)*VEC2(1)
            VEC3(3)= VEC1(1)*VEC2(2)-VEC1(2)*VEC2(1)
            C1=(XYZ(3*(NEWATOM-1)+1)-XYZ(3*(CONLIST(1)-1)+1))*VEC1(1)+ &
  &            (XYZ(3*(NEWATOM-1)+2)-XYZ(3*(CONLIST(1)-1)+2))*VEC1(2)+ &
  &            (XYZ(3*(NEWATOM-1)+3)-XYZ(3*(CONLIST(1)-1)+3))*VEC1(3)
            C2=(XYZ(3*(NEWATOM-1)+1)-XYZ(3*(CONLIST(1)-1)+1))*VEC2(1)+ &
  &            (XYZ(3*(NEWATOM-1)+2)-XYZ(3*(CONLIST(1)-1)+2))*VEC2(2)+ &
  &            (XYZ(3*(NEWATOM-1)+3)-XYZ(3*(CONLIST(1)-1)+3))*VEC2(3)
            C3=(XYZ(3*(NEWATOM-1)+1)-XYZ(3*(CONLIST(1)-1)+1))*VEC3(1)+ &
  &            (XYZ(3*(NEWATOM-1)+2)-XYZ(3*(CONLIST(1)-1)+2))*VEC3(2)+ &
  &            (XYZ(3*(NEWATOM-1)+3)-XYZ(3*(CONLIST(1)-1)+3))*VEC3(3)
            DO J1=2,INTIMAGE+1
               VEC1(1:3)=XYZ((J1-1)*3*NATOMS+3*(CONLIST(2)-1)+1:(J1-1)*3*NATOMS+3*(CONLIST(2)-1)+3) &
  &                     -XYZ((J1-1)*3*NATOMS+3*(CONLIST(1)-1)+1:(J1-1)*3*NATOMS+3*(CONLIST(1)-1)+3)
               DUMMY=SQRT(VEC1(1)**2+VEC1(2)**2+VEC1(3)**2)
               IF (DUMMY.NE.0.0D0) VEC1(1:3)=VEC1(1:3)/DUMMY
               VEC2(1:3)=XYZ((J1-1)*3*NATOMS+3*(CONLIST(3)-1)+1:(J1-1)*3*NATOMS+3*(CONLIST(3)-1)+3) &
  &                     -XYZ((J1-1)*3*NATOMS+3*(CONLIST(1)-1)+1:(J1-1)*3*NATOMS+3*(CONLIST(1)-1)+3)
               DUMMY=VEC1(1)*VEC2(1)+VEC1(2)*VEC2(2)+VEC1(3)*VEC2(3)
               VEC2(1:3)=VEC2(1:3)-DUMMY*VEC1(1:3)
               DUMMY=SQRT(VEC2(1)**2+VEC2(2)**2+VEC2(3)**2)
               IF (DUMMY.NE.0.0D0) VEC2(1:3)=VEC2(1:3)/DUMMY
               VEC3(1)= VEC1(2)*VEC2(3)-VEC1(3)*VEC2(2)
               VEC3(2)=-VEC1(1)*VEC2(3)+VEC1(3)*VEC2(1)
               VEC3(3)= VEC1(1)*VEC2(2)-VEC1(2)*VEC2(1)
               XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)= &
  &            XYZ((J1-1)*3*NATOMS+3*(CONLIST(1)-1)+1:(J1-1)*3*NATOMS+3*(CONLIST(1)-1)+3)+C1*VEC1(1:3)+C2*VEC2(1:3)+C3*VEC3(1:3)
            ENDDO
            CALL CHECKREP(INTIMAGE,XYZ,NOPT) ! set up repulsive neighbour list
            IF (CHECKCONINT) THEN
               CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
            ELSE
               CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
            ENDIF
            ESAVE0=ETOTAL
            DO J1=2,INTIMAGE+1
               XSAVE0(1:3,J1)=XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)
            ENDDO
         ENDIF
         IF (NDFORNEWATOM.GE.3) THEN
!
! Choose three atoms from the BESTPRESERVEDN list at random with bias towards the 
! start of the list. Let the relative weight for position i be 1/i**2 and calculate
! the sum to normalise.
!
            DUMMY=0.0D0
            DO J1=1,NDFORNEWATOM
!              DUMMY=DUMMY+1.0D0/(1.0D0*J1)
!              DUMMY=DUMMY+1.0D0/(1.0D0*BESTPRESERVEDD(J1))
               DUMMY=DUMMY+1.0D0/(1.0D0*J1**2)
            ENDDO
            N1=0; N2=0; N3=0
            DO WHILE (N3.EQ.0)
               DUMMY2=0.0D0
               RAN1=DPRAND()*DUMMY
               DO J1=1,NDFORNEWATOM
!                 DUMMY2=DUMMY2+1.0D0/(1.0D0*J1)
!                 DUMMY2=DUMMY2+1.0D0/(1.0D0*BESTPRESERVEDD(J1))
                  DUMMY2=DUMMY2+1.0D0/(1.0D0*J1**2)
                  IF (DUMMY2.GE.RAN1) THEN
                     IF ((J1.EQ.N1).OR.(J1.EQ.N2)) EXIT ! already chosen
                     IF (N1.EQ.0) THEN
                        N1=J1
                        EXIT
                     ENDIF
                     IF (N2.EQ.0) THEN
                        N2=J1
                        EXIT
                     ENDIF
                     N3=J1
                     EXIT
                  ENDIF
               ENDDO
            ENDDO
            IF (DEBUG) PRINT '(A,3I6,A)',' intlbfgs> choosing positions ',N1,N2,N3,' in best preserved list'
            IF (DEBUG) PRINT '(A,3I6)',' intlbfgs> atoms are ',BESTPRESERVEDN(N1),BESTPRESERVEDN(N2),BESTPRESERVEDN(N3)
!           IF (DEBUG) PRINT '(A,3I6,A)',' intlbfgs> full list has length ',NDFORNEWATOM
!           IF (DEBUG) PRINT '(20I6)',BESTPRESERVEDN(1:NDFORNEWATOM)

!
! Move the new atom consistently in the local environment of the three active atoms with the
! best preserved absolute distances or the shortest average distances in the end points.
! Check the energies and compare linear interpolation as well, then choose the interpolation
! with the lowest energy.
! Make a local orthogonal coordinate system and use constant components in this basis.
!
            VEC1(1:3)=XYZ(3*(BESTPRESERVEDN(N2)-1)+1:3*(BESTPRESERVEDN(N2)-1)+3) &
  &                  -XYZ(3*(BESTPRESERVEDN(N1)-1)+1:3*(BESTPRESERVEDN(N1)-1)+3)
            DUMMY=SQRT(VEC1(1)**2+VEC1(2)**2+VEC1(3)**2)
            IF (DUMMY.NE.0.0D0) VEC1(1:3)=VEC1(1:3)/DUMMY
            VEC2(1:3)=XYZ(3*(BESTPRESERVEDN(N3)-1)+1:3*(BESTPRESERVEDN(N3)-1)+3) &
  &                  -XYZ(3*(BESTPRESERVEDN(N1)-1)+1:3*(BESTPRESERVEDN(N1)-1)+3)
            DUMMY=VEC1(1)*VEC2(1)+VEC1(2)*VEC2(2)+VEC1(3)*VEC2(3)
            VEC2(1:3)=VEC2(1:3)-DUMMY*VEC1(1:3)
            DUMMY=SQRT(VEC2(1)**2+VEC2(2)**2+VEC2(3)**2)
            IF (DUMMY.NE.0.0D0) VEC2(1:3)=VEC2(1:3)/DUMMY
            VEC3(1)= VEC1(2)*VEC2(3)-VEC1(3)*VEC2(2)
            VEC3(2)=-VEC1(1)*VEC2(3)+VEC1(3)*VEC2(1)
            VEC3(3)= VEC1(1)*VEC2(2)-VEC1(2)*VEC2(1)
            C1=(XYZ(3*(NEWATOM-1)+1)-XYZ(3*(BESTPRESERVEDN(N1)-1)+1))*VEC1(1)+ &
  &            (XYZ(3*(NEWATOM-1)+2)-XYZ(3*(BESTPRESERVEDN(N1)-1)+2))*VEC1(2)+ &
  &            (XYZ(3*(NEWATOM-1)+3)-XYZ(3*(BESTPRESERVEDN(N1)-1)+3))*VEC1(3)
            C2=(XYZ(3*(NEWATOM-1)+1)-XYZ(3*(BESTPRESERVEDN(N1)-1)+1))*VEC2(1)+ &
  &            (XYZ(3*(NEWATOM-1)+2)-XYZ(3*(BESTPRESERVEDN(N1)-1)+2))*VEC2(2)+ &
  &            (XYZ(3*(NEWATOM-1)+3)-XYZ(3*(BESTPRESERVEDN(N1)-1)+3))*VEC2(3)
            C3=(XYZ(3*(NEWATOM-1)+1)-XYZ(3*(BESTPRESERVEDN(N1)-1)+1))*VEC3(1)+ &
  &            (XYZ(3*(NEWATOM-1)+2)-XYZ(3*(BESTPRESERVEDN(N1)-1)+2))*VEC3(2)+ &
  &            (XYZ(3*(NEWATOM-1)+3)-XYZ(3*(BESTPRESERVEDN(N1)-1)+3))*VEC3(3)
            DO J1=2,INTIMAGE+1
               VEC1(1:3)=XYZ((J1-1)*3*NATOMS+3*(BESTPRESERVEDN(N2)-1)+1:(J1-1)*3*NATOMS+3*(BESTPRESERVEDN(N2)-1)+3) &
  &                     -XYZ((J1-1)*3*NATOMS+3*(BESTPRESERVEDN(N1)-1)+1:(J1-1)*3*NATOMS+3*(BESTPRESERVEDN(N1)-1)+3)
               DUMMY=SQRT(VEC1(1)**2+VEC1(2)**2+VEC1(3)**2)
               IF (DUMMY.NE.0.0D0) VEC1(1:3)=VEC1(1:3)/DUMMY
               VEC2(1:3)=XYZ((J1-1)*3*NATOMS+3*(BESTPRESERVEDN(N3)-1)+1:(J1-1)*3*NATOMS+3*(BESTPRESERVEDN(N3)-1)+3) &
  &                     -XYZ((J1-1)*3*NATOMS+3*(BESTPRESERVEDN(N1)-1)+1:(J1-1)*3*NATOMS+3*(BESTPRESERVEDN(N1)-1)+3)
               DUMMY=VEC1(1)*VEC2(1)+VEC1(2)*VEC2(2)+VEC1(3)*VEC2(3)
               VEC2(1:3)=VEC2(1:3)-DUMMY*VEC1(1:3)
               DUMMY=SQRT(VEC2(1)**2+VEC2(2)**2+VEC2(3)**2)
               IF (DUMMY.NE.0.0D0) VEC2(1:3)=VEC2(1:3)/DUMMY
               VEC3(1)= VEC1(2)*VEC2(3)-VEC1(3)*VEC2(2)
               VEC3(2)=-VEC1(1)*VEC2(3)+VEC1(3)*VEC2(1)
               VEC3(3)= VEC1(1)*VEC2(2)-VEC1(2)*VEC2(1)
               XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)= &
  &            XYZ((J1-1)*3*NATOMS+3*(BESTPRESERVEDN(N1)-1)+1:(J1-1)*3*NATOMS+3*(BESTPRESERVEDN(N1)-1)+3)+ &
  &                   C1*VEC1(1:3)+C2*VEC2(1:3)+C3*VEC3(1:3)
            ENDDO

            CALL CHECKREP(INTIMAGE,XYZ,NOPT) ! set up repulsive neighbour list
            IF (CHECKCONINT) THEN
               CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
            ELSE
               CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
            ENDIF
            ESAVED=ETOTAL
            DO J1=2,INTIMAGE+1
               XSAVED(1:3,J1)=XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)
            ENDDO
         ENDIF

         IF (NCFORNEWATOM.GE.3) THEN
!
! Choose three atoms from the BESTCLOSEST list at random with bias towards the
! start of the list. Let the relative weight for position i be 1/i**2 and calculate
! the sum to normalise.
!
            DUMMY=0.0D0
            DO J1=1,NCFORNEWATOM
!              DUMMY=DUMMY+1.0D0/(1.0D0*J1)
!              DUMMY=DUMMY+1.0D0/(1.0D0*BESTCLOSESTD(J1))
               DUMMY=DUMMY+1.0D0/(1.0D0*J1**2)
            ENDDO
            N1=0; N2=0; N3=0
            DO WHILE (N3.EQ.0)
               DUMMY2=0.0D0
               RAN1=DPRAND()*DUMMY
               DO J1=1,NCFORNEWATOM
!                 DUMMY2=DUMMY2+1.0D0/(1.0D0*J1)
!                 DUMMY2=DUMMY2+1.0D0/(1.0D0*BESTCLOSESTD(J1))
                  DUMMY2=DUMMY2+1.0D0/(1.0D0*J1**2)
                  IF (DUMMY2.GE.RAN1) THEN
                     IF ((J1.EQ.N1).OR.(J1.EQ.N2)) EXIT ! already chosen
                     IF (N1.EQ.0) THEN
                        N1=J1
                        EXIT
                     ENDIF
                     IF (N2.EQ.0) THEN
                        N2=J1
                        EXIT
                     ENDIF
                     N3=J1
                     EXIT
                  ENDIF
               ENDDO
            ENDDO
            IF (DEBUG) PRINT '(A,3I6,A)',' intlbfgs> choosing positions ',N1,N2,N3,' in closest list'

            VEC1(1:3)=XYZ(3*(BESTCLOSESTN(N2)-1)+1:3*(BESTCLOSESTN(N2)-1)+3)-XYZ(3*(BESTCLOSESTN(N1)-1)+1:3*(BESTCLOSESTN(N1)-1)+3)
            DUMMY=SQRT(VEC1(1)**2+VEC1(2)**2+VEC1(3)**2)
            IF (DUMMY.NE.0.0D0) VEC1(1:3)=VEC1(1:3)/DUMMY
            VEC2(1:3)=XYZ(3*(BESTCLOSESTN(N3)-1)+1:3*(BESTCLOSESTN(N3)-1)+3)-XYZ(3*(BESTCLOSESTN(N1)-1)+1:3*(BESTCLOSESTN(N1)-1)+3)
            DUMMY=VEC1(1)*VEC2(1)+VEC1(2)*VEC2(2)+VEC1(3)*VEC2(3)
            VEC2(1:3)=VEC2(1:3)-DUMMY*VEC1(1:3)
            DUMMY=SQRT(VEC2(1)**2+VEC2(2)**2+VEC2(3)**2)
            IF (DUMMY.NE.0.0D0) VEC2(1:3)=VEC2(1:3)/DUMMY
            VEC3(1)= VEC1(2)*VEC2(3)-VEC1(3)*VEC2(2)
            VEC3(2)=-VEC1(1)*VEC2(3)+VEC1(3)*VEC2(1)
            VEC3(3)= VEC1(1)*VEC2(2)-VEC1(2)*VEC2(1)
            C1=(XYZ(3*(NEWATOM-1)+1)-XYZ(3*(BESTCLOSESTN(N1)-1)+1))*VEC1(1)+ &
  &            (XYZ(3*(NEWATOM-1)+2)-XYZ(3*(BESTCLOSESTN(N1)-1)+2))*VEC1(2)+ &
  &            (XYZ(3*(NEWATOM-1)+3)-XYZ(3*(BESTCLOSESTN(N1)-1)+3))*VEC1(3)
            C2=(XYZ(3*(NEWATOM-1)+1)-XYZ(3*(BESTCLOSESTN(N1)-1)+1))*VEC2(1)+ &
  &            (XYZ(3*(NEWATOM-1)+2)-XYZ(3*(BESTCLOSESTN(N1)-1)+2))*VEC2(2)+ &
  &            (XYZ(3*(NEWATOM-1)+3)-XYZ(3*(BESTCLOSESTN(N1)-1)+3))*VEC2(3)
            C3=(XYZ(3*(NEWATOM-1)+1)-XYZ(3*(BESTCLOSESTN(N1)-1)+1))*VEC3(1)+ &
  &            (XYZ(3*(NEWATOM-1)+2)-XYZ(3*(BESTCLOSESTN(N1)-1)+2))*VEC3(2)+ &
  &            (XYZ(3*(NEWATOM-1)+3)-XYZ(3*(BESTCLOSESTN(N1)-1)+3))*VEC3(3)
            DO J1=2,INTIMAGE+1
               VEC1(1:3)=XYZ((J1-1)*3*NATOMS+3*(BESTCLOSESTN(N2)-1)+1:(J1-1)*3*NATOMS+3*(BESTCLOSESTN(N2)-1)+3) &
  &                     -XYZ((J1-1)*3*NATOMS+3*(BESTCLOSESTN(N1)-1)+1:(J1-1)*3*NATOMS+3*(BESTCLOSESTN(N1)-1)+3)
               DUMMY=SQRT(VEC1(1)**2+VEC1(2)**2+VEC1(3)**2)
               IF (DUMMY.NE.0.0D0) VEC1(1:3)=VEC1(1:3)/DUMMY
               VEC2(1:3)=XYZ((J1-1)*3*NATOMS+3*(BESTCLOSESTN(N3)-1)+1:(J1-1)*3*NATOMS+3*(BESTCLOSESTN(N3)-1)+3) &
  &                     -XYZ((J1-1)*3*NATOMS+3*(BESTCLOSESTN(N1)-1)+1:(J1-1)*3*NATOMS+3*(BESTCLOSESTN(N1)-1)+3)
               DUMMY=VEC1(1)*VEC2(1)+VEC1(2)*VEC2(2)+VEC1(3)*VEC2(3)
               VEC2(1:3)=VEC2(1:3)-DUMMY*VEC1(1:3)
               DUMMY=SQRT(VEC2(1)**2+VEC2(2)**2+VEC2(3)**2)
               IF (DUMMY.NE.0.0D0) VEC2(1:3)=VEC2(1:3)/DUMMY
               VEC3(1)= VEC1(2)*VEC2(3)-VEC1(3)*VEC2(2)
               VEC3(2)=-VEC1(1)*VEC2(3)+VEC1(3)*VEC2(1)
               VEC3(3)= VEC1(1)*VEC2(2)-VEC1(2)*VEC2(1)
               XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)= &
  &            XYZ((J1-1)*3*NATOMS+3*(BESTCLOSESTN(N1)-1)+1:(J1-1)*3*NATOMS+3*(BESTCLOSESTN(N1)-1)+3)+ &
  &                   C1*VEC1(1:3)+C2*VEC2(1:3)+C3*VEC3(1:3)
            ENDDO

            CALL CHECKREP(INTIMAGE,XYZ,NOPT) ! set up repulsive neighbour list
            IF (CHECKCONINT) THEN
               CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
            ELSE
               CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
            ENDIF
            ESAVEC=ETOTAL
            DO J1=2,INTIMAGE+1
               XSAVEC(1:3,J1)=XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)
            ENDDO
         ENDIF
!
! Standard linear interpolation, with constraint distance scaled by FRAC.
! Works for FRAC as small as 0.1 with repulsion turned off.
! We use an appropriately weighted displacement from atom CONLIST(1) using the displacements
! in the two end points.
!
         FRAC=1.0D0
         DO J1=2,INTIMAGE+1
            XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1)=XYZ((J1-1)*3*NATOMS+3*(CONLIST(1)-1)+1)  &
 &            +(INTIMAGE-J1+2)*FRAC*(XYZ(3*(NEWATOM-1)+1)-XYZ(3*(CONLIST(1)-1)+1))/(INTIMAGE+1) &
 &   +(J1-1)*(XYZ(3*NATOMS*(INTIMAGE+1)+3*(NEWATOM-1)+1)-XYZ(3*NATOMS*(INTIMAGE+1)+3*(CONLIST(1)-1)+1))/(INTIMAGE+1)
            XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+2)=XYZ((J1-1)*3*NATOMS+3*(CONLIST(1)-1)+2)  &
 &            +(INTIMAGE-J1+2)*FRAC*(XYZ(3*(NEWATOM-1)+2)-XYZ(3*(CONLIST(1)-1)+2))/(INTIMAGE+1) &
 &   +(J1-1)*(XYZ(3*NATOMS*(INTIMAGE+1)+3*(NEWATOM-1)+2)-XYZ(3*NATOMS*(INTIMAGE+1)+3*(CONLIST(1)-1)+2))/(INTIMAGE+1)
            XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+3)=XYZ((J1-1)*3*NATOMS+3*(CONLIST(1)-1)+3)  &
 &            +(INTIMAGE-J1+2)*FRAC*(XYZ(3*(NEWATOM-1)+3)-XYZ(3*(CONLIST(1)-1)+3))/(INTIMAGE+1) &
 &   +(J1-1)*(XYZ(3*NATOMS*(INTIMAGE+1)+3*(NEWATOM-1)+3)-XYZ(3*NATOMS*(INTIMAGE+1)+3*(CONLIST(1)-1)+3))/(INTIMAGE+1)
         ENDDO
         CALL CHECKREP(INTIMAGE,XYZ,NOPT) ! set up repulsive neighbour list
         IF (CHECKCONINT) THEN
            CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
         ELSE
            CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
         ENDIF
         IF (DEBUG) PRINT '(A,4G15.5)',' intlbfgs> energies for constrained, preserved, closest, and linear schemes=', &
  &                 ESAVE0,ESAVED,ESAVEC,ETOTAL
         IF ((ETOTAL.LT.ESAVEC).AND.(ETOTAL.LT.ESAVED).AND.(ETOTAL.LT.ESAVE0)) THEN
            IF (DEBUG) PRINT '(A,2G20.10)',' intlbfgs> lowest energy from linear interpolation'
         ELSE IF ((ESAVEC.LT.ESAVED).AND.(ESAVEC.LT.ESAVE0)) THEN
            IF (DEBUG) PRINT '(A,2G20.10)',' intlbfgs> lowest energy from interpolation using closest atoms'
            DO J1=2,INTIMAGE+1
               XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)=XSAVEC(1:3,J1)
            ENDDO
            ETOTAL=ESAVEC
         ELSE IF (ESAVED.LT.ESAVE0) THEN
            IF (DEBUG) PRINT '(A,2G20.10)',' intlbfgs> lowest energy from interpolation using preserved distances'
            DO J1=2,INTIMAGE+1
               XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)=XSAVED(1:3,J1)
            ENDDO
            ETOTAL=ESAVED
         ELSE 
            IF (DEBUG) PRINT '(A,2G20.10)',' intlbfgs> lowest energy from interpolation using closest constraints'
            DO J1=2,INTIMAGE+1
               XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)=XSAVE0(1:3,J1)
            ENDDO
            ETOTAL=ESAVE0
         ENDIF
      ENDIF
      NADDED=NADDED+1
      IF (NADDED.LT.NTOADD) GOTO 542
      CALL CHECKREP(INTIMAGE,XYZ,NOPT) ! set up repulsive neighbour list
!
! need a new gradient since the active atom has changed !
!
      IF (CHECKCONINT) THEN
         CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
      ELSE
         CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
      ENDIF
      NLASTGOODE=NITERDONE
      LASTGOODE=ETOTAL
      GOODSAVE(1:D)=X(1:D)
   ENDIF 

   MAIN: IF (NITERDONE==1) THEN
        POINT = 0
        DIAG(1:D)=INTDGUESS
        SEARCHSTEP(0,1:D)= -G(1:D)*DIAG                           ! NR STEP FOR DIAGONAL INVERSE HESSIAN
        GTMP(1:D)        = SEARCHSTEP(0,1:D)
        GNORM            = MAX(SQRT(DOT_PRODUCT(G(1:D),G(1:D))),1.0D-100)
        STP(1:D)         = MIN(1.0D0/GNORM, GNORM)           ! MAKE THE FIRST GUESS FOR THE STEP LENGTH CAUTIOUS
   ELSE MAIN
        BOUND=NITERDONE-1
        IF (NITERDONE.GT.INTMUPDATE) BOUND=INTMUPDATE
        YS=DOT_PRODUCT( GDIF(NPT/D,:), SEARCHSTEP(NPT/D,:)  )
        IF (YS==0.0D0) YS=1.0D0
    
        ! Update estimate of diagonal inverse Hessian elements.
        ! We divide by both YS and YY at different points, so they had better not be zero!

        YY=DOT_PRODUCT( GDIF(NPT/D,:) , GDIF(NPT/D,:) )
        IF (YY==0.0D0) YY=1.0D0
!       DIAG = ABS(YS/YY)
        DIAG = YS/YY
      
        ! COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980, "Updating quasi-Newton matrices with limited storage",
        ! Mathematics of Computation, Vol.35, No.151, pp. 773-782
        CP= POINT; IF (POINT==0) CP = INTMUPDATE
        RHO1(CP)=1.0D0/YS
        GTMP(1:D) = -G(1:D)
        CP= POINT 
                   
        DO I= 1,BOUND 
!            CP = CP - 1; IF (CP == -1) CP = M - 1
             CP = CP - 1; IF (CP == -1) CP = INTMUPDATE - 1
             SQ= DOT_PRODUCT( SEARCHSTEP(CP,1:D),GTMP(1:D) )
             ALPHA(CP+1) = RHO1(CP+1) * SQ
             GTMP(1:D)        = -ALPHA(CP+1)*GDIF(CP,1:D) + GTMP(1:D)
        ENDDO
              
        GTMP(1:D)=DIAG*GTMP(1:D)

        DO I=1,BOUND
             YR= DOT_PRODUCT( GDIF(CP,1:D) , GTMP )
             BETA= RHO1(CP+1)*YR
             BETA= ALPHA(CP+1)-BETA
             GTMP(1:D) = BETA*SEARCHSTEP(CP,1:D) + GTMP(1:D)
             CP=CP+1
!            IF (CP==M) CP=0
             IF (CP==INTMUPDATE) CP=0
        ENDDO
              
        STP(1:D) = 1.0D0
   ENDIF MAIN

   !  Store the new search direction
   IF (NITERDONE.GT.1) SEARCHSTEP(POINT,1:D)=GTMP(1:D)
      
!
! If the number of images has changed since G was declared then G is not the same
! size as Gtmp and Dot_Product cannot be used.
!
!  IF (Dot_Product(G,Gtmp)/SQRT( Dot_Product(G,G)*Dot_Product(Gtmp,Gtmp) ) > 0.0D0) THEN
!
!  Separate sqrt;s to avoid overflow.
!
   IF (DDOT(D,G,1,GTMP,1)/MAX(1.0D-100,SQRT( DDOT(D,G,1,G,1))*SQRT(DDOT(D,GTMP,1,GTMP,1)) ) > 0.0D0) THEN
        IF (DEBUG) PRINT*,'Search direction has positive projection onto gradient - reversing step'
        GTMP(1:D)=-GTMP(1:D)
        SEARCHSTEP(POINT,1:D)=GTMP(1:D)
   ENDIF
   GTMP(1:D)=G(1:D)

!  We should apply the maximum LBFGS step to each image separately.
!  However, using different scale factors for different images leads to huge
!  discontinuities! Now take the minimum scale factor for all images. DJW 26/11/07

   STPMIN=1.0D0
   DO J2=1,INTIMAGE
      STEPIMAGE(J2) = SQRT(DOT_PRODUCT(SEARCHSTEP(POINT,NOPT*(J2-1)+1:NOPT*J2),SEARCHSTEP(POINT,NOPT*(J2-1)+1:NOPT*J2)))
      DUMMY=STEPIMAGE(J2)
      IF (STEPIMAGE(J2) > MAXINTBFGS) THEN
           STP(NOPT*(J2-1)+1:NOPT*J2) = MAXINTBFGS/STEPIMAGE(J2)
           STPMIN=MIN(STPMIN,STP(NOPT*(J2-1)+1))
      ENDIF
!     PRINT '(A,I8,3G20.10)',' image,initial step size,STP,prod=',J2,DUMMY,STP(NOPT*(J2-1)+1),STEPIMAGE(J2)*STP(NOPT*(J2-1)+1)
   ENDDO
   STP(1:D)=STPMIN

! EFK: decide whether to freeze some nodes
   IF (FREEZENODEST) THEN
      TOTGNORM=SQRT(DOT_PRODUCT(G(1:NOPT*INTIMAGE),G(1:NOPT*INTIMAGE))/INTIMAGE)
      NIMAGEFREEZE=0
      DO IM=1,INTIMAGE
         TESTG=SQRT(DOT_PRODUCT(G(NOPT*(IM-1)+1:NOPT*IM),G(NOPT*(IM-1)+1:NOPT*IM)))
         IMGFREEZE(IM)=.FALSE.
         IF (TOTGNORM.NE.0.0D0) THEN
            IF (TESTG/TOTGNORM.LT.FREEZETOL) THEN
!              IF (DEBUG) PRINT '(A,I6,2G20.10)', ' intlbfgs> Freezing image: ', IM, TESTG, TOTGNORM
               IMGFREEZE(IM)=.TRUE.
               STEPIMAGE(IM)=0.0D0
               NIMAGEFREEZE=NIMAGEFREEZE+1
               STP(NOPT*(IM-1)+1:NOPT*IM)=0.0D0
            ENDIF
         ENDIF
      ENDDO
      IF (DEBUG) PRINT '(2(A,I6))', ' intlbfgs> Number of frozen images=',NIMAGEFREEZE,' / ',INTIMAGE
   ENDIF
   !  We now have the proposed step - update geometry and calculate new gradient
   NDECREASE=0
20 X(1:D) = X(1:D) + STP(1:D)*SEARCHSTEP(POINT,1:D)

   IF (.NOT.SWITCHED) THEN
!     IF ((RMS.LT.INTRMSTOL*1.0D10).AND.(MOD(NITERDONE,10).EQ.0).AND.(NSTEPSMAX-NITERDONE.GT.100)) &
! &               CALL CHECKSEP(NMAXINT,NMININT,INTIMAGE,XYZ,NOPT,NATOMS)
      IF (MOD(NITERDONE,50).EQ.0) THEN
         DMAX=0.0D0
         DMIN=HUGE(1.0D0)
         DO J1=1,INTIMAGE+1
            DUMMY=0.0D0
            DO J2=1,3*NATOMS
               IF (ATOMACTIVE((J2-1)/3+1)) THEN
                  DUMMY=DUMMY+( XYZ((J1-1)*3*NATOMS+J2) - XYZ(J1*3*NATOMS+J2) )**2
               ENDIF
            ENDDO
            DUMMY=SQRT(DUMMY)
            IF (DUMMY.GT.DMAX) THEN
               DMAX=DUMMY
               JMAX=J1
            ENDIF
            IF (DUMMY.LT.DMIN) THEN
               DMIN=DUMMY
               JMIN=J1
            ENDIF
            IF (DEBUG) PRINT '(A,I6,A,I6,A,G20.10)',' intlbfgs> distance between images ', &
  &                                                  J1,' and ',J1+1,' is ',DUMMY
         ENDDO
         IF ((DMAX.GT.IMSEPMAX).AND.(INTIMAGE.LT.MAXINTIMAGE)) THEN
            PRINT '(A,I6,A,I6)',' intlbfgs> Add an image between ',JMAX,' and ',JMAX+1
            ALLOCATE(DPTMP(3*NATOMS*(INTIMAGE+2)))
            DPTMP(1:3*NATOMS*(INTIMAGE+2))=XYZ(1:3*NATOMS*(INTIMAGE+2))
            DEALLOCATE(XYZ)
            ALLOCATE(XYZ(3*NATOMS*(INTIMAGE+3)))
            XYZ(1:3*NATOMS*JMAX)=DPTMP(1:3*NATOMS*JMAX)
            XYZ(3*NATOMS*JMAX+1:3*NATOMS*(JMAX+1))=(DPTMP(3*NATOMS*(JMAX-1)+1:3*NATOMS*JMAX) &
  &                                               + DPTMP(3*NATOMS*JMAX+1:3*NATOMS*(JMAX+1)))/2.0D0
            XYZ(3*NATOMS*(JMAX+1)+1:3*NATOMS*(INTIMAGE+3))=DPTMP(3*NATOMS*JMAX+1:3*NATOMS*(INTIMAGE+2))
!
! Save step-taking memories in SEARCHSTEP and GDIF.
! These arrays run from 0 to INTMUPDATE over memories and
! 1:NOPT*INTIMAGE over only the variable images.
!
            DEALLOCATE(DPTMP)
            ALLOCATE(D2TMP(0:INTMUPDATE,1:NOPT*INTIMAGE))
            D2TMP(0:INTMUPDATE,1:NOPT*INTIMAGE)=SEARCHSTEP(0:INTMUPDATE,1:NOPT*INTIMAGE)
            DEALLOCATE(SEARCHSTEP)
            ALLOCATE(SEARCHSTEP(0:INTMUPDATE,1:NOPT*(INTIMAGE+1)))
            DO J1=1,INTMUPDATE
               IF (JMAX.GT.1) SEARCHSTEP(J1,1:3*NATOMS*(JMAX-1))=D2TMP(J1,1:3*NATOMS*(JMAX-1))
               IF (JMAX.LT.INTIMAGE+1) SEARCHSTEP(J1,3*NATOMS*JMAX+1:3*NATOMS*(INTIMAGE+1))= &
  &                 D2TMP(J1,3*NATOMS*(JMAX-1)+1:3*NATOMS*INTIMAGE)
               SEARCHSTEP(J1,3*NATOMS*(JMAX-1)+1:3*NATOMS*JMAX)= &
  &                             D2TMP(J1,3*NATOMS*(MIN(JMAX,INTIMAGE)-1)+1:3*NATOMS*MIN(JMAX,INTIMAGE))
            ENDDO
            D2TMP(0:INTMUPDATE,1:NOPT*INTIMAGE)=GDIF(0:INTMUPDATE,1:NOPT*INTIMAGE)
            DEALLOCATE(GDIF)
            ALLOCATE(GDIF(0:INTMUPDATE,1:NOPT*(INTIMAGE+1)))
            DO J1=1,INTMUPDATE
               IF (JMAX.GT.1) GDIF(J1,1:3*NATOMS*(JMAX-1))=D2TMP(J1,1:3*NATOMS*(JMAX-1))
               IF (JMAX.LT.INTIMAGE+1) GDIF(J1,3*NATOMS*JMAX+1:3*NATOMS*(INTIMAGE+1))= &
  &                 D2TMP(J1,3*NATOMS*(JMAX-1)+1:3*NATOMS*INTIMAGE)
               GDIF(J1,3*NATOMS*(JMAX-1)+1:3*NATOMS*JMAX)= &
  &                       D2TMP(J1,3*NATOMS*(MIN(JMAX,INTIMAGE)-1)+1:3*NATOMS*MIN(JMAX,INTIMAGE))
            ENDDO
            DEALLOCATE(D2TMP)

            DEALLOCATE(SAVEX,TRUEEE,GOODSAVE,EEETMP,MYGTMP,XSAVED,XSAVEC,XSAVE0,GTMP,GGG, &
  &                    DIAG,STP,GLAST,XSAVE,EEE,STEPIMAGE,CHECKG,IMGFREEZE)
            ALLOCATE(SAVEX(3*NATOMS*(INTIMAGE+1)),TRUEEE(INTIMAGE+3),GOODSAVE(3*NATOMS*(INTIMAGE+3)), &
  &                  EEETMP(INTIMAGE+3), MYGTMP(3*NATOMS*(INTIMAGE+1)), XSAVED(3,INTIMAGE+3), &
  &                  XSAVEC(3,INTIMAGE+3), XSAVE0(3,INTIMAGE+3), GTMP(3*NATOMS*(INTIMAGE+1)), &
  &                  DIAG(3*NATOMS*(INTIMAGE+1)), STP(3*NATOMS*(INTIMAGE+1)), &
  &                  GLAST(NOPT*(INTIMAGE+1)), &
  &                  XSAVE(NOPT*(INTIMAGE+1)), CHECKG(NOPT*(INTIMAGE+1)), IMGFREEZE(INTIMAGE+1), &
  &                  EEE(INTIMAGE+3), STEPIMAGE(INTIMAGE+1), GGG(3*NATOMS*(INTIMAGE+3)))
            GGG(1:3*NATOMS*(INTIMAGE+3))=0.0D0
            SAVEX(1:3*NATOMS*(INTIMAGE+1))=0.0D0
            TRUEEE(1:INTIMAGE+3)=0.0D0
            EEETMP(1:INTIMAGE+3)=0.0D0
            MYGTMP(1:3*NATOMS*(INTIMAGE+1))=0.0D0
            XSAVED(1:3,1:INTIMAGE+3)=0.0D0
            XSAVEC(1:3,1:INTIMAGE+3)=0.0D0
            XSAVE0(1:3,1:INTIMAGE+3)=0.0D0
            GTMP(1:3*NATOMS*(INTIMAGE+1))=0.0D0
            DIAG(1:3*NATOMS*(INTIMAGE+1))=0.0D0
            STP(1:3*NATOMS*(INTIMAGE+1))=0.0D0
            GLAST(1:NOPT*(INTIMAGE+1))=0.0D0
            XSAVE(1:NOPT*(INTIMAGE+1))=0.0D0
            CHECKG(1:NOPT*(INTIMAGE+1))=.FALSE.
            IMGFREEZE(1:INTIMAGE+1)=.FALSE.
            EEE(1:INTIMAGE+3)=0.0D0
            STEPIMAGE(1:INTIMAGE+1)=0.0D0

            X=>XYZ(NOPT+1:NOPT*(INTIMAGE+2))
            G=>GGG(NOPT+1:NOPT*(INTIMAGE+2))
            INTIMAGE=INTIMAGE+1
            D=NOPT*INTIMAGE
            CALL CHECKREP(INTIMAGE,XYZ,NOPT)
            IF (CHECKCONINT) THEN
               CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
            ELSE
               CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
            ENDIF
         ELSEIF ((DMIN.LT.IMSEPMIN).AND.(INTIMAGE.GT.1)) THEN
            IF (JMIN.EQ.1) JMIN=2
            PRINT '(A,I6,A,I6)',' intlbfgs> Remove image ',JMIN
            ALLOCATE(DPTMP(3*NATOMS*(INTIMAGE+2)))
            DPTMP(1:3*NATOMS*(INTIMAGE+2))=XYZ(1:3*NATOMS*(INTIMAGE+2))
            DEALLOCATE(XYZ)
            ALLOCATE(XYZ(3*NATOMS*(INTIMAGE+1)))
            XYZ(1:3*NATOMS*(JMIN-1))=DPTMP(1:3*NATOMS*(JMIN-1))
            XYZ(3*NATOMS*(JMIN-1)+1:3*NATOMS*(INTIMAGE+1))=DPTMP(3*NATOMS*JMIN+1:3*NATOMS*(INTIMAGE+2))

            DEALLOCATE(DPTMP)
!
! Save step-taking memories in SEARCHSTEP and GDIF.
! These arrays run from 0 to INTMUPDATE over memories and
! 1:NOPT*INTIMAGE over only the variable images.
!
            ALLOCATE(D2TMP(0:INTMUPDATE,1:NOPT*INTIMAGE))
            D2TMP(0:INTMUPDATE,1:NOPT*INTIMAGE)=SEARCHSTEP(0:INTMUPDATE,1:NOPT*INTIMAGE)
            DEALLOCATE(SEARCHSTEP)
            ALLOCATE(SEARCHSTEP(0:INTMUPDATE,1:NOPT*(INTIMAGE-1)))
            DO J1=1,INTMUPDATE
               SEARCHSTEP(J1,1:3*NATOMS*(JMIN-2))=D2TMP(J1,1:3*NATOMS*(JMIN-2))
               SEARCHSTEP(J1,3*NATOMS*(JMIN-2)+1:3*NATOMS*(INTIMAGE-1))= &
  &                     D2TMP(J1,3*NATOMS*(JMIN-1)+1:3*NATOMS*INTIMAGE)
            ENDDO
            D2TMP(0:INTMUPDATE,1:NOPT*INTIMAGE)=GDIF(0:INTMUPDATE,1:NOPT*INTIMAGE)
            DEALLOCATE(GDIF)
            ALLOCATE(GDIF(0:INTMUPDATE,1:NOPT*(INTIMAGE-1)))
            DO J1=1,INTMUPDATE
               GDIF(J1,1:3*NATOMS*(JMIN-2))=D2TMP(J1,1:3*NATOMS*(JMIN-2))
               GDIF(J1,3*NATOMS*(JMIN-2)+1:3*NATOMS*(INTIMAGE-1))= &
  &                     D2TMP(J1,3*NATOMS*(JMIN-1)+1:3*NATOMS*INTIMAGE)
            ENDDO
            DEALLOCATE(D2TMP)

            DEALLOCATE(SAVEX,TRUEEE,GOODSAVE,EEETMP,MYGTMP,XSAVED,XSAVEC,XSAVE0,GTMP,GGG, &
  &                    DIAG,STP,GLAST,XSAVE,EEE,STEPIMAGE,CHECKG,IMGFREEZE)
            ALLOCATE(SAVEX(3*NATOMS*(INTIMAGE-1)),TRUEEE(INTIMAGE+1),GOODSAVE(3*NATOMS*(INTIMAGE+1)), &
  &                  EEETMP(INTIMAGE+1), MYGTMP(3*NATOMS*(INTIMAGE-1)), XSAVED(3,INTIMAGE+1), &
  &                  XSAVEC(3,INTIMAGE+1), XSAVE0(3,INTIMAGE+1), GTMP(3*NATOMS*(INTIMAGE-1)), &
  &                  DIAG(3*NATOMS*(INTIMAGE-1)), STP(3*NATOMS*(INTIMAGE-1)), &
  &                  GLAST(NOPT*(INTIMAGE-1)), &
  &                  XSAVE(NOPT*(INTIMAGE-1)), CHECKG(NOPT*(INTIMAGE-1)), IMGFREEZE(INTIMAGE-1), &
  &                  EEE(INTIMAGE+1), STEPIMAGE(INTIMAGE-1), GGG(3*NATOMS*(INTIMAGE+1)))
            GGG(1:3*NATOMS*(INTIMAGE+1))=0.0D0
            SAVEX(1:3*NATOMS*(INTIMAGE-1))=0.0D0
            TRUEEE(1:INTIMAGE+1)=0.0D0
            EEETMP(1:INTIMAGE+1)=0.0D0
            MYGTMP(1:3*NATOMS*(INTIMAGE-1))=0.0D0
            XSAVED(1:3,1:INTIMAGE+1)=0.0D0
            XSAVEC(1:3,1:INTIMAGE+1)=0.0D0
            XSAVE0(1:3,1:INTIMAGE+1)=0.0D0
            GTMP(1:3*NATOMS*(INTIMAGE-1))=0.0D0
            DIAG(1:3*NATOMS*(INTIMAGE-1))=0.0D0
            STP(1:3*NATOMS*(INTIMAGE-1))=0.0D0
            GLAST(1:NOPT*(INTIMAGE-1))=0.0D0
            XSAVE(1:NOPT*(INTIMAGE-1))=0.0D0
            CHECKG(1:NOPT*(INTIMAGE-1))=.FALSE.
            IMGFREEZE(1:INTIMAGE-1)=.FALSE.
            EEE(1:INTIMAGE+1)=0.0D0
            STEPIMAGE(1:INTIMAGE-1)=0.0D0

            X=>XYZ(NOPT+1:NOPT*(INTIMAGE))
            G=>GGG(NOPT+1:NOPT*(INTIMAGE))
            INTIMAGE=INTIMAGE-1
            D=NOPT*INTIMAGE
            CALL CHECKREP(INTIMAGE,XYZ,NOPT)
            IF (CHECKCONINT) THEN
               CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
            ELSE
               CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
            ENDIF
         ENDIF
      ELSE
         IF (MOD(NITERDONE,25).EQ.0) CALL CHECKREP(INTIMAGE,XYZ,NOPT)
         IF (CHECKCONINT) THEN
            CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
         ELSE
            CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
         ENDIF
      ENDIF
      IF ((ETOTAL-EOLD.LT.1.0D100).OR.ADDATOM) THEN ! MAXERISE effectively set to 1.0D100 here
!     IF ((ETOTAL-EOLD.LT.2.0D-1).OR.ADDATOM) THEN 
         EOLD=ETOTAL
         GLAST(1:D)=G(1:D)
         XSAVE(1:D)=X(1:D)
      ELSE
         NDECREASE=NDECREASE+1
         IF (NDECREASE.GT.5) THEN
            NFAIL=NFAIL+1
            WRITE(*,'(A,I6)') ' intlbfgs> WARNING *** in lbfgs cannot find a lower energy, NFAIL=',NFAIL
            X(1:D)=XSAVE(1:D)
            G(1:D)=GLAST(1:D)
         ELSE
            X(1:D)=XSAVE(1:D)
            G(1:D)=GLAST(1:D)
            STP(1:D)=STP(1:D)/10.0D0
            WRITE(*,'(A,G25.15,A,G25.15,A)') ' intlbfgs> energy increased from ',EOLD,' to ',ETOTAL, &
     &          ' decreasing step size'
            GOTO 20
         ENDIF
      ENDIF
      ADDATOM=.FALSE.
   ELSE ! combine constraint and true potentials
!     IF ((RMS.LT.INTRMSTOL*1.0D10).AND.(MOD(NITERDONE,10).EQ.0).AND.(NSTEPSMAX-NITERDONE.GT.100)) &
! &               CALL CHECKSEP(NMAXINT,NMININT,INTIMAGE,XYZ,NOPT)

!!!
!
! Check that MAKE_CONPOT produces the same constraints and repulsions - this is to debug MAKE_CONPOT
!
!     MINCOORDS(1,1:NOPT)=XYZ(1:NOPT)
!     MINCOORDS(2,1:NOPT)=XYZ(NOPT*(INTIMAGE+1)+1:NOPT*(INTIMAGE+2))
!     PRINT '(A)',' intlbfgs> Before make_conpot'
!     CALL CHECKREP(INTIMAGE,XYZ,NOPT)
!     DO J2=1,NCONSTRAINT
!        PRINT '(A,I6,L5,2I6,2F20.10)','J2,CONACTIVE,CONI,CONJ,CONDISTREF,CONDISTREFLOCAL=', &
! &                      J2,CONACTIVE(J2),CONI(J2),CONJ(J2),CONDISTREF(J2),CONDISTREFLOCAL(J2)
!     ENDDO
!     DO J2=1,NREPULSIVE
!        PRINT '(A,3I6,F20.10)','J2,REPI,REPJ,REPCUT=',J2,REPI(J2),REPJ(J2),REPCUT(J2)
!     ENDDO
!     DO J2=1,NNREPULSIVE
!        PRINT '(A,3I6,F20.10)','J2,NREPI,NREPJ,NREPCUT=',J2,NREPI(J2),NREPJ(J2),NREPCUT(J2)
!     ENDDO
!     PRINT '(A)',' intlbfgs> Calling make_conpot'
!     CALL MAKE_CONPOT(2,MINCOORDS)
!     CALL CHECKREP(INTIMAGE,XYZ,NOPT)
!     DO J2=1,NCONSTRAINT
!        PRINT '(A,I6,L5,2I6,2F20.10)','J2,CONACTIVE,CONI,CONJ,CONDISTREF,CONDISTREFLOCAL=', &
! &                      J2,CONACTIVE(J2),CONI(J2),CONJ(J2),CONDISTREF(J2),CONDISTREFLOCAL(J2)
!     ENDDO
!     DO J2=1,NREPULSIVE
!        PRINT '(A,3I6,F20.10)','J2,REPI,REPJ,REPCUT=',J2,REPI(J2),REPJ(J2),REPCUT(J2)
!     ENDDO
!     DO J2=1,NNREPULSIVE
!        PRINT '(A,3I6,F20.10)','J2,NREPI,NREPJ,NREPCUT=',J2,NREPI(J2),NREPJ(J2),NREPCUT(J2)
!     ENDDO
!     STOP
!!! DJW
      IF (MOD(NITERDONE,25).EQ.0) CALL CHECKREP(INTIMAGE,XYZ,NOPT)
      ETOTALTMP=0.0D0
      DO J4=2,INTIMAGE+1
         CALL POTENTIAL(XYZ(NOPT*(J4-1)+1:NOPT*J4),EEE(J4),GGG(NOPT*(J4-1)+1:NOPT*J4),.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
         ETOTALTMP=ETOTALTMP+EEE(J4)
      ENDDO
      RMSTMP=RMS
      EEETMP(1:INTIMAGE+2)=EEE(1:INTIMAGE+2)
      MYGTMP(1:D)=G(1:D)
      IF (CHECKCONINT) THEN
         CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
      ELSE
         CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
      ENDIF
      ETOTAL=USEFRAC*ETOTALTMP+(1.0D0-USEFRAC)*ETOTAL
      RMS=USEFRAC*RMSTMP+(1.0D0-USEFRAC)*RMS
      G(1:D)=USEFRAC*MYGTMP(1:D)+(1.0D0-USEFRAC)*G(1:D)
      EEE(1:INTIMAGE+2)=USEFRAC*EEETMP(1:INTIMAGE+2)+(1.0D0-USEFRAC)*EEE(1:INTIMAGE+2)
!     USEFRAC=USEFRAC+INTCONFRAC
      IF (USEFRAC.GE.1.0D0) PRINT '(A,I6)',' intlbfgs> switching off constraint potential completely at step ',NITERDONE
   ENDIF
   IF (ETOTAL/INTIMAGE.LT.COLDFUSIONLIMIT) THEN
      WRITE(*,'(A,2G20.10)') ' intlbfgs> Cold fusion diagnosed - step discarded, energy, limit=',ETOTAL/INTIMAGE,COLDFUSIONLIMIT
      DEALLOCATE(CONI,CONJ,CONDISTREF,REPI,REPJ,NREPI,NREPJ,REPCUT,NREPCUT)
      DEALLOCATE(SAVEX,TRUEEE,GOODSAVE, EEETMP, MYGTMP, XSAVED, XSAVEC, XSAVE0, GTMP, &
  &              DIAG, STP, SEARCHSTEP, GDIF,GLAST, XSAVE, XYZ, GGG, CHECKG, IMGFREEZE, EEE, STEPIMAGE)
      INTIMAGE=INTIMAGESAVE
      NTSFOUND=0
      NMINFOUND=0
      RETURN
   ENDIF

   STEPTOT = SUM(STEPIMAGE)/INTIMAGE

   IF (DEBUG) THEN
      WRITE(*,'(A,I6,2G20.10,G20.10,I8)') ' intlbfgs> steps: ',NITERDONE,ETOTAL/INTIMAGE,RMS,STEPTOT,NACTIVE
      CALL FLUSH(6,ISTAT)
   ENDIF

   IF (.NOT.SWITCHED) THEN
      IF ((NITERDONE-NLASTGOODE.GT.INTRELSTEPS).AND.((ETOTAL.GT.LASTGOODE).OR.(ETOTAL/INTIMAGE.GT.MAXCONE*1.0D8))) THEN
         PRINT '(2(A,I6))',' intlbfgs> Backtracking ',NBACKTRACK,' steps, active atoms=',NACTIVE-NBACKTRACK
         NTRIES(NEWATOM)=NTRIES(NEWATOM)+1
         IF (FREEZENODEST) IMGFREEZE(1:INTIMAGE)=.FALSE.
!
! Backtrack by removing the last NBACKTRACK atoms along with their active constraints and
! repulsions.
!
         DO J1=1,NBACKTRACK
            NDUMMY=TURNONORDER(NACTIVE-J1+1)
            IF (DEBUG) PRINT '(A,I6,A,2I6)',' intlbfgs> Turning off active atom ',NDUMMY
            DO J2=1,NCONSTRAINT
               IF (.NOT.CONACTIVE(J2)) CYCLE
               IF ((CONI(J2).EQ.NDUMMY).OR.(CONJ(J2).EQ.NDUMMY)) THEN
                  CONACTIVE(J2)=.FALSE.
                  IF (DEBUG) PRINT '(A,I6,A,2I6)',' intlbfgs> Turning off constraint ',J2,' for atoms ',CONI(J2),CONJ(J2)
               ENDIF
            ENDDO
            ATOMACTIVE(NDUMMY)=.FALSE.
         ENDDO
!
! Reconstruct repulsions. 
!
         NREPULSIVE=0
         DO J1=1,NATOMS
            IF (.NOT.ATOMACTIVE(J1)) CYCLE ! identify active atoms
            DO J2=J1+1,NATOMS
               IF (.NOT.ATOMACTIVE(J2)) CYCLE ! identify active atoms
               IF (ABS(J1-J2).LE.INTREPSEP) CYCLE ! no repulsion for atoms too close in sequence
               DO J3=1,NCONSTRAINT
                  IF (.NOT.CONACTIVE(J3)) CYCLE ! identify active constraints 
                  IF (((CONI(J3).EQ.J1).AND.(CONJ(J3).EQ.J2)).OR. &
  &                   ((CONJ(J3).EQ.J1).AND.(CONI(J3).EQ.J2))) GOTO 548
               ENDDO
               DMIN=1.0D100
               DO J3=1,INTIMAGE+2,INTIMAGE+1 ! only consider the end-point distances
                  DF=SQRT((XYZ((J3-1)*3*NATOMS+3*(J2-1)+1)-XYZ((J3-1)*3*NATOMS+3*(J1-1)+1))**2+ &
  &                       (XYZ((J3-1)*3*NATOMS+3*(J2-1)+2)-XYZ((J3-1)*3*NATOMS+3*(J1-1)+2))**2+ &
  &                       (XYZ((J3-1)*3*NATOMS+3*(J2-1)+3)-XYZ((J3-1)*3*NATOMS+3*(J1-1)+3))**2)
                  IF (DF.LT.DMIN) DMIN=DF
               ENDDO
!
! Use the minimum of the end point distances and INTCONSTRAINREPCUT for each contact.
!
               DMIN=MIN(DMIN-1.0D-3,INTCONSTRAINREPCUT)
               NREPULSIVE=NREPULSIVE+1
               REPI(NREPULSIVE)=J1
               REPJ(NREPULSIVE)=J2
               REPCUT(NREPULSIVE)=DMIN
548            CONTINUE
            ENDDO
         ENDDO

         NACTIVE=NACTIVE-NBACKTRACK
         NBACKTRACK=MIN(MIN(1.0D0*(NBACKTRACK+1),1.0D0*20),0.1D0*(NACTIVE-2))
         IF (DEBUG) PRINT '(A,I6)',' intlbfgs> Number of atoms to backtrack is now ',NBACKTRACK
         NDUMMY=0
         DO J1=1,NATOMS
            IF (ATOMACTIVE(J1)) NDUMMY=NDUMMY+1
         ENDDO
         IF (NDUMMY.NE.NACTIVE) THEN
            PRINT '(A,I6)',' intlbfgs> ERROR *** inconsistency in number of active atoms. Should be ',NACTIVE
            DO J1=1,NATOMS
               IF (ATOMACTIVE(J1)) PRINT '(A,I6)',' active atom ',J1
            ENDDO
            STOP
         ENDIF
         ADDATOM=.TRUE.

         CALL CHECKREP(INTIMAGE,XYZ,NOPT)
         IF (CHECKCONINT) THEN
            CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
         ELSE
            CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
         ENDIF
      ENDIF
      LASTGOODE=ETOTAL
   ENDIF

   EXITSTATUS=0
   IF (INTIMAGE.EQ.0) THEN ! sanity check - should be unnecessary
      PRINT '(A,2I10)','intlbfgs> ERROR *** INTIMAGE,SIZE(DIAG,1)=',INTIMAGE,SIZE(DIAG,1)
      STOP
   ENDIF
   INTDGUESS=DIAG(1) ! should be ok for subsequent runs of the same system DJW
   IF (RMS<=INTRMSTOL.AND.NITERDONE>1) EXITSTATUS=1
   IF (NITERDONE==NSTEPSMAX) EXITSTATUS=2

   IF (.FALSE.) THEN
      CHECKG(1:D)=.FALSE.
      DO J1=1,D
         IF (ABS(G(J1)).GT.1.0D-6) THEN
            PRINT '(3I6,G20.10)',J1,2+(J1-1)/(3*NATOMS),(J1-3*NATOMS*((J1-1)/(3*NATOMS))-1)/3+1,G(J1)
            CHECKG(J1)=.TRUE.
         ENDIF
     ENDDO
!!!!!!!!!!!!!!!!!!!
!     NDUMMY=NREPULSIVE
!     NCONSTRAINT=0
!     NREPULSIVE=0
!!!!!!!!!!!!!!!!!!!
      IF (CHECKCONINT) THEN
         CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
      ELSE
         CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
      ENDIF
      GLAST(1:D)=G(1:D)
      DIFF=1.0D-6
      PRINT '(A,I6)',' intlbfgs> analytic and numerical gradients: D=',D
      DO J2=1,D
         IF (.NOT.CHECKG(J2)) CYCLE
         X(J2)=X(J2)+DIFF
!        PRINT '(A,I6)',' intlbfgs> calling congrad + for coordinate J2'
         IF (CHECKCONINT) THEN
            CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
         ELSE
            CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
         ENDIF
         EPLUS=ETOTAL
         X(J2)=X(J2)-2.0D0*DIFF
!        PRINT '(A,I6)',' intlbfgs> calling congrad - for coordinate J2'
         IF (CHECKCONINT) THEN
            CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
         ELSE
            CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
         ENDIF
         EMINUS=ETOTAL
         X(J2)=X(J2)+DIFF
         IF (ABS(GLAST(J2)).NE.0.0D0) THEN
            IF (100.0D0*ABS((GLAST(J2)-(EPLUS-EMINUS)/(2.0D0*DIFF))/GLAST(J2)).GT.10.0D0) THEN
               WRITE(*,'(A,3I8,3G20.10)') 'error ',(J2-1)/NOPT+1,(J2-NOPT*((J2-1)/NOPT)-1)/3+1,J2, &
  &                                 GLAST(J2),(EPLUS-EMINUS)/(2.0D0*DIFF), &
  &                                 (EPLUS-EMINUS)/(2.0D0*DIFF*GLAST(J2))
            ELSE
               WRITE(*,'(A,3I8,3G20.10)') 'OK    ',(J2-1)/NOPT+1,(J2-NOPT*((J2-1)/NOPT)-1)/3+1,J2, &
  &                                       GLAST(J2),(EPLUS-EMINUS)/(2.0D0*DIFF), &
  &                                       (EPLUS-EMINUS)/(2.0D0*DIFF*GLAST(J2))
            ENDIF
         ENDIF
      ENDDO
   ENDIF

   IF (EXITSTATUS > 0) THEN  
      IF ((.NOT.SWITCHED).AND.(EXITSTATUS.EQ.1)) THEN ! add active atom or restart with true potential on
         IF (ETOTAL/INTIMAGE.GT.MAXCONE) GOTO 777
         IF (NACTIVE.LT.NATOMS) THEN 
            ADDATOM=.TRUE.
            GOTO 777
         ENDIF
         CALL MYCPU_TIME(FTIME,.FALSE.)
         PRINT '(A,I6,A,F12.6,A,I6,A,F10.1)',' intlbfgs> switch on true potential at step ',NITERDONE, &
  &                                     ' fraction=',INTCONFRAC,' images=',INTIMAGE,' time=',FTIME-STIME
         PRINT '(A,I6,A,F15.6)',' intlbfgs> Allowing ',INTCONSTEPS,' further optimization steps'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   IF (.FALSE.) THEN
      CHECKG(1:D)=.FALSE.
      DO J1=1,D
         IF (ABS(G(J1)).GT.1.0D-6) THEN
            PRINT '(3I6,G20.10)',J1,2+(J1-1)/(3*NATOMS),(J1-3*NATOMS*((J1-1)/(3*NATOMS))-1)/3+1,G(J1)
            CHECKG(J1)=.TRUE.
         ENDIF
     ENDDO
      IF (CHECKCONINT) THEN
         CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
      ELSE
         CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
      ENDIF
      GLAST(1:D)=G(1:D)
      DIFF=1.0D-6
      PRINT '(A,I6)',' intlbfgs> analytic and numerical gradients: D=',D
      DO J2=1,D
         IF (.NOT.CHECKG(J2)) CYCLE
         X(J2)=X(J2)+DIFF
!        PRINT '(A,I6)',' intlbfgs> calling congrad + for coordinate J2'
         IF (CHECKCONINT) THEN
            CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
         ELSE
            CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
         ENDIF
         EPLUS=ETOTAL
         X(J2)=X(J2)-2.0D0*DIFF
!        PRINT '(A,I6)',' intlbfgs> calling congrad - for coordinate J2'
         IF (CHECKCONINT) THEN
            CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
         ELSE
            CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
         ENDIF
         EMINUS=ETOTAL
         X(J2)=X(J2)+DIFF
         IF (ABS(GLAST(J2)).NE.0.0D0) THEN
            IF (100.0D0*ABS((GLAST(J2)-(EPLUS-EMINUS)/(2.0D0*DIFF))/GLAST(J2)).GT.10.0D0) THEN
               WRITE(*,'(A,3I8,3G20.10)') 'error ',(J2-1)/NOPT+1,(J2-NOPT*((J2-1)/NOPT)-1)/3+1,J2, &
  &                                 GLAST(J2),(EPLUS-EMINUS)/(2.0D0*DIFF), &
  &                                 (EPLUS-EMINUS)/(2.0D0*DIFF*GLAST(J2))
            ELSE
               WRITE(*,'(A,3I8,3G20.10)') 'OK    ',(J2-1)/NOPT+1,(J2-NOPT*((J2-1)/NOPT)-1)/3+1,J2, &
  &                                       GLAST(J2),(EPLUS-EMINUS)/(2.0D0*DIFF), &
  &                                       (EPLUS-EMINUS)/(2.0D0*DIFF*GLAST(J2))
            ENDIF
         ENDIF
      ENDDO
   ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         DO J1=1,NATOMS
            IF (.NOT.ATOMACTIVE(J1)) THEN
               PRINT '(A,I6,A,I6,A)',' intlbfgs> ERROR *** number of active atoms=',NACTIVE,' but atom ',J1,' is not active'
            ENDIF
         ENDDO
         NSTEPSMAX=NITERDONE+INTCONSTEPS
         SWITCHED=.TRUE.
         USEFRAC=INTCONFRAC
         GOTO 777
      ELSEIF ((.NOT.SWITCHED).AND.(EXITSTATUS.EQ.2)) THEN 
         PRINT '(A,I6)',' intlbfgs> ERROR *** number of active atoms at final step=',NACTIVE
         STOP
      ELSEIF (DEBUG) THEN
         PRINT '(A,I6,A,I6)','intlbfgs> energies for images:'
         PRINT '(I6,F20.10)',(J2,EEE(J2),J2=1,INTIMAGE+2)
      ENDIF
      EXIT
   ENDIF
   777 CONTINUE
!
! Compute the new step and gradient change
!
   NPT=POINT*D
   SEARCHSTEP(POINT,:) = STP*SEARCHSTEP(POINT,:)
   GDIF(POINT,:)=G-GTMP
   POINT=POINT+1; IF (POINT==INTMUPDATE) POINT=0

   IF (DUMPINTXYZ.AND.MOD(NITERDONE,DUMPINTXYZFREQ)==0) CALL RWG(NITERDONE,INTIMAGE,XYZ)
   IF (DUMPINTEOS.AND.MOD(NITERDONE,DUMPINTEOSFREQ)==0) CALL WRITEPROFILE(NITERDONE,EEE,INTIMAGE)
   PREVGRAD=RMS

   NITERDONE=NITERDONE+1
   IF (NITERDONE.GT.NSTEPSMAX) EXIT

ENDDO ! end of main do loop over counter NITERDONE

IF (.NOT.SWITCHED) THEN 
   PRINT '(A,I6)',' intlbfgs> ERROR *** number of active atoms at final step=',NACTIVE,' no potential switch'
   STOP
ENDIF
IF (EXITSTATUS.EQ.1) THEN
   WRITE(*,'(A,I6,A,G20.10,A,G20.10)') ' intlbfgs> Converged after ',NITERDONE,' steps, energy per image=',ETOTAL/INTIMAGE, &
  &                               ' RMS gradient=',RMS
ELSEIF (EXITSTATUS.EQ.2) THEN
   WRITE(*,'(A,I6,A,G20.10,A,G20.10)') ' intlbfgs> After ',NITERDONE,' steps, energy per image=',ETOTAL/INTIMAGE, &
  &                               ' RMS gradient=',RMS
ENDIF
!
! Linear interpolation for constraint potential and real potential separately.
! Constraint potential need not be flat if we have done some steps with both
! potentials turned on.
!
DINCREMENT=0.02D0
DTOTAL=0.0D0
OPEN(UNIT=753,FILE='intenergy',STATUS='UNKNOWN')
!
! local maxima must have NSIDE higher energies on each side
! This has the desirable side-effect that we don't bother with
! images that are essentially collapsed on each other - their
! spacing will probably be < DINCREMENT, or 5*DINCREMENT.
!
NSIDE=5
INTTST=.TRUE. ! try passing local maxima back as ts guesses
! INTTST=.FALSE. 
NTSFOUND=0
NMINFOUND=0
PRINTOPTIMIZETS=DEBUG
DO J1=1,INTIMAGE+1
   DUMMY=0.0D0
   DO J2=1,3*NATOMS
      DUMMY=DUMMY+( XYZ((J1-1)*3*NATOMS+J2) - XYZ(J1*3*NATOMS+J2) )**2
   ENDDO
   DUMMY=SQRT(DUMMY)
   DIST=0.0D0
   PRINT '(A,I6,A,I6,A,G20.10)',' intlbfgs> distance between images ',J1,' and ',J1+1,' is ',DUMMY
   NDUMMY=DUMMY/DINCREMENT+1
   ALLOCATE(EINT(NDUMMY))
   J3=1

   intloop: DO
      LOCALCOORDS(1:3*NATOMS)=((DUMMY-DIST)*XYZ((J1-1)*3*NATOMS+1:J1*3*NATOMS)+ &
  &                                    DIST*XYZ(J1*3*NATOMS+1:(J1+1)*3*NATOMS))/DUMMY
      CALL POTENTIAL(LOCALCOORDS,EREAL,VNEW,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
      If (DEBUG) PRINT '(A,3G20.10)',' intlbfgs> ',DTOTAL+DIST,EREAL
      WRITE(753,'(3G20.10)') DTOTAL+DIST,EREAL
      DIST=DIST+DINCREMENT
      EINT(J3)=EREAL
      IF (INTTST) THEN
         IF (J3-NSIDE.GT.0) THEN
            DO J4=MAX(J3-2*NSIDE,1),J3
               IF (J4.EQ.J3-NSIDE) CYCLE
               IF (EINT(J3-NSIDE).LT.EINT(J4)) GOTO 432
            ENDDO
!
! We have a ts candidate. Try optimising it!
!
            CALL MYCPU_TIME(STARTTIME,.FALSE.)
            KNOWG=.FALSE.
            KNOWE=.FALSE. ! to be safe!
            LOCALCOORDS(1:NOPT)= &
  &                 ((DUMMY-(J3-NSIDE-1)*DINCREMENT)*XYZ((J1-1)*NOPT+1:J1*NOPT)+ &
  &                         (J3-NSIDE-1)*DINCREMENT *XYZ(J1*NOPT+1:(J1+1)*NOPT))/DUMMY
            IF (BFGSTST) THEN
               VECS(1:NOPT)=(XYZ((J1-1)*NOPT+1:J1*NOPT)-XYZ(J1*NOPT+1:(J1+1)*NOPT))/DUMMY
               CALL BFGSTS(NSTEPS,LOCALCOORDS,  &
  &               EDUMMY,LGDUMMY,TSCONVERGED,RMS,EVALMIN,EVALMAX,VECS,ITDONE,.TRUE.,PRINTOPTIMIZETS)
            ELSE
               CALL EFOL(LOCALCOORDS,TSCONVERGED,NSTEPS,EDUMMY,ITDONE,EVALMIN,DEBUG,XDIAG,2)
            ENDIF
            CALL MYCPU_TIME(TIME0,.FALSE.)
            IF (TSCONVERGED) THEN
               NTSFOUND=NTSFOUND+1
!
! Save coordinates and direction vector between images to use as starting guess
! for the eigenvector.
!
               ALLOCATE(TSFOUND(NTSFOUND)%E,TSFOUND(NTSFOUND)%COORD(NOPT), &
  &                     TSFOUND(NTSFOUND)%EVALMIN,TSFOUND(NTSFOUND)%VECS(NOPT))
               TSFOUND(NTSFOUND)%VECS(1:NOPT)=VECS(1:NOPT)
               TSFOUND(NTSFOUND)%COORD(1:NOPT)=LOCALCOORDS(1:NOPT)
               TSFOUND(NTSFOUND)%E=EDUMMY
               TSFOUND(NTSFOUND)%EVALMIN=EVALMIN
               PRINT '(A,I6,A,G20.10,A,F10.1)',' intlbfgs> transition state found, iterations=',ITDONE, &
  &                                  ' energy=',EDUMMY,' time=',TIME0-STARTTIME
            ENDIF
432         CONTINUE
         ENDIF
      ENDIF
      J3=J3+1
      IF (DIST.GT.DUMMY) EXIT INTLOOP
      IF (J3.GT.NDUMMY) THEN
         PRINT '(A,I6)',' intlbfgs> ERROR *** number of interpolated energies should not be ',J3
      ENDIF
   ENDDO intloop
   DTOTAL=DTOTAL+DUMMY
   DEALLOCATE(EINT)
ENDDO

LOCALCOORDS(1:3*NATOMS)=XYZ((INTIMAGE+1)*3*NATOMS+1:(INTIMAGE+2)*3*NATOMS)
CALL POTENTIAL(LOCALCOORDS,EREAL,VNEW,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
PRINT '(A,3G20.10)',' intlbfgs> ',DTOTAL,EREAL
WRITE(753,'(3G20.10)') DTOTAL,EREAL
CLOSE(753)

IF (.NOT.INTTST) THEN
   PTEST=.FALSE.
   INTTST=.FALSE.
   PRINT '(A,I8)',' intlbfgs> minimising all the images - results written to images.min'
   OPEN(987,FILE='images.min',STATUS='UNKNOWN')
   WRITE(987,'(I6)') NATOMS
   WRITE(987,'(A)') 'start - image 1'
   WRITE(987,'(A,3G20.10)') (ZSYM(J2),XYZ(3*(J2-1)+1:3*(J2-1)+3),J2=1,NATOMS)
   DO J1=2,INTIMAGE+1
      KNOWG=.FALSE.
      KNOWE=.FALSE. ! could use EEE value
      PRINT '(A,I8,A,F20.10)',' intlbfgs> minimising image ',J1,' initial energy=',EEE(J1)

!     BSMIN=.TRUE.
!     DEBUG=.TRUE.
!     CALL ODESD(100,XYZ(NOPT*(J1-1)+1:NOPT*J1),MFLAG,ITDONE,.TRUE.)
!     DEBUG=.FALSE.
!     BSMIN=.FALSE.

!     KNOWG=.FALSE.
!     KNOWE=.FALSE. 
!     PTEST=.TRUE.
      CALL MYLBFGS(NOPT,MUPDATE,XYZ(NOPT*(J1-1)+1:NOPT*J1),.FALSE., &
   &               MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,.TRUE.,ITDONE,PTEST,VNEW,.TRUE.,.FALSE.)
!     PTEST=.FALSE.
   
      IF (MFLAG) THEN
         NMINFOUND=NMINFOUND+1
!
!  We have to communicate the minima found back to tryconnect using the data structure
!  set up for new transition states.
!  Added new variable MINFOUND to allow for this check in tryconnect.
!  It seems impossible to make intlbfgs see isnewmin and addnewmin for some reason.
!
         ALLOCATE(MINFOUND(NMINFOUND)%E,MINFOUND(NMINFOUND)%COORD(NOPT))
         MINFOUND(NMINFOUND)%COORD(1:NOPT)=XYZ(NOPT*(J1-1)+1:NOPT*J1)
         MINFOUND(NMINFOUND)%E=EREAL
         WRITE(987,'(I6)') NATOMS
         WRITE(987,'(A,I5)') 'image ',J1
         WRITE(987,'(A,3G20.10)') (ZSYM(J2), MINFOUND(NMINFOUND)%COORD(3*(J2-1)+1:3*(J2-1)+3),J2=1,NATOMS)
      ENDIF
   ENDDO
   WRITE(987,'(I6)') NATOMS
   WRITE(987,'(A)') 'finish - image INTIMAGE+2'
   WRITE(987,'(A,3G20.10)') (ZSYM(J2),XYZ(NOPT*(INTIMAGE+1)+3*(J2-1)+1:NOPT*(INTIMAGE+1)+3*(J2-1)+3),J2=1,NATOMS)
   CLOSE(987)
ENDIF
 
DEALLOCATE(CONI,CONJ,CONDISTREF,REPI,REPJ,NREPI,NREPJ,REPCUT,NREPCUT)
DEALLOCATE(SAVEX,TRUEEE,GOODSAVE, EEETMP, MYGTMP, XSAVED, XSAVEC, XSAVE0, GTMP, &
  &      DIAG, STP, SEARCHSTEP, GDIF,GLAST, XSAVE, XYZ, GGG, CHECKG, IMGFREEZE, EEE, STEPIMAGE)
INTIMAGE=INTIMAGESAVE

END SUBROUTINE INTLBFGS
!
! Possible redistribution of images for INTCONSTRAINT depending upon distances.
!
SUBROUTINE CHECKSEP(NMAXINT,NMININT,INTIMAGE,XYZ,NOPT,NATOMS)
IMPLICIT NONE
INTEGER NSEPMAX, NSEPMIN, J, NMININT, NMAXINT, INTIMAGE, NOPT, J1, J2, NATOMS
DOUBLE PRECISION SEPMAX, SEPMIN, XYZ(*), DUMMY

RETURN !!! DJW

IF ((NMININT.EQ.NMAXINT).OR.(NMININT.EQ.NMAXINT+1)) THEN
   PRINT '(A,2I6)',' checksep> skipping image redistribution for images ',NMININT,NMAXINT
   RETURN
ENDIF
IF ((NMININT.EQ.1).OR.(NMININT.EQ.INTIMAGE+2)) THEN
   PRINT '(A,I6)',' checksep> ERROR *** NMININT=',NMININT
ENDIF
IF ((NMAXINT.EQ.1).OR.(NMAXINT.EQ.INTIMAGE+2)) THEN
   PRINT '(A,I6)',' checksep> ERROR *** NMAXINT=',NMAXINT
ENDIF
! 
! DVEC(J) contains the distance between image J and image J+1
!
!      SEPMAX=-1.0D0
!      SEPMIN=1.0D100
!      DO J=1,INTIMAGE+1
!         IF (DVEC(J).GT.SEPMAX) THEN
!            SEPMAX=DVEC(J)
!            NSEPMAX=J
!         ENDIF
!      ENDDO
!      DO J=2,INTIMAGE+1
!         IF (DVEC(J-1)+DVEC(J).LT.SEPMIN) THEN
!            SEPMIN=DVEC(J-1)+DVEC(J)
!            NSEPMIN=J
!         ENDIF
!      ENDDO
!      PRINT '(A,F20.10,A,I6,A,I6)',' checksep> maximum image separation=',SEPMAX,' for images ',NSEPMAX,' and ',NSEPMAX+1
!      PRINT '(A,F20.10,A,I6)',' checksep> minimum sum of image separations=',SEPMIN,' for image ',NSEPMIN

!    IF (SEPMIN*2.0D0.LT.SEPMAX) THEN ! redistribute images

IF (.TRUE.) THEN ! redistribute images
   PRINT '(A,I6,A,2I6)',' checksep> removing image ',NMININT,' and adding one between images ',NMAXINT,NMAXINT+1
!  IF (NSEPMIN.LT.NSEPMAX) THEN
   IF (NMININT.LT.NMAXINT) THEN
!     DO J=NSEPMIN,NSEPMAX-1 ! move image J+1 to position J for images J=NSEPMIN+1 to NSEPMAX-1
      DO J=NMININT,NMAXINT-1 ! move image J+1 to position J for images J=NMININT+1 to NMAXINT-1
         XYZ(NOPT*(J-1)+1:NOPT*J)=XYZ(NOPT*J+1:NOPT*(J+1))
      ENDDO
!     XYZ(NOPT*(NSEPMAX-1)+1:NOPT*NSEPMAX)=(XYZ(NOPT*(NSEPMAX-1)+1:NOPT*NSEPMAX)+XYZ(NOPT*NSEPMAX+1:NOPT*(NSEPMAX+1)))/2.0D0
      XYZ(NOPT*(NMAXINT-1)+1:NOPT*NMAXINT)=(XYZ(NOPT*(NMAXINT-1)+1:NOPT*NMAXINT)+XYZ(NOPT*NMAXINT+1:NOPT*(NMAXINT+1)))/2.0D0
   ELSE
!     DO J=NSEPMIN,NSEPMAX+2,-1 ! move image J-1 to position J for images J=NSEPMIN-1 to NSEPMAX+1
      DO J=NMININT,NMAXINT+2,-1 ! move image J-1 to position J for images J=NMININT-1 to NMAXINT+1
         PRINT '(2(A,I6))',' putting image ',J-1,' in image ',J
         XYZ(NOPT*(J-1)+1:NOPT*J)=XYZ(NOPT*(J-2)+1:NOPT*(J-1))
      ENDDO
      XYZ(NOPT*NMAXINT+1:NOPT*(NMAXINT+1))=(XYZ(NOPT*NMAXINT+1:NOPT*(NMAXINT+1))+XYZ(NOPT*(NMAXINT-1)+1:NOPT*NMAXINT))/2.0D0
   ENDIF
ENDIF

END SUBROUTINE CHECKSEP
!
! Neighbour list for repulsions to reduce cost of constraint potential.
!
SUBROUTINE CHECKREP(INTIMAGE,XYZ,NOPT)
USE KEY,ONLY : NREPI, NREPJ, NREPCUT, NNREPULSIVE, NREPULSIVE, REPI, REPJ, REPCUT, DEBUG
USE PORFUNCS
IMPLICIT NONE
INTEGER JJ, KK, NI1, NJ1, NI2, NJ2, INTIMAGE, NOPT, ISTAT
DOUBLE PRECISION LDIST, XYZ(NOPT*(INTIMAGE+2))
DOUBLE PRECISION R1AX,R1AY,R1AZ,R2AX,R2AY,R2AZ,R1BX,R1BY,R1BZ,R2BX,R2BY,R2BZ,DMIN
LOGICAL NOINT
 
NNREPULSIVE=0
DO JJ=1,NREPULSIVE
 CALL FLUSH(6,ISTAT)
   DO KK=1,INTIMAGE+2 ! first check for standard distances within threshold
      LDIST=SQRT((XYZ((KK-1)*NOPT+3*(REPI(JJ)-1)+1)-XYZ((KK-1)*NOPT+3*(REPJ(JJ)-1)+1))**2 &
  &             +(XYZ((KK-1)*NOPT+3*(REPI(JJ)-1)+2)-XYZ((KK-1)*NOPT+3*(REPJ(JJ)-1)+2))**2 &
  &             +(XYZ((KK-1)*NOPT+3*(REPI(JJ)-1)+3)-XYZ((KK-1)*NOPT+3*(REPJ(JJ)-1)+3))**2)
      IF (LDIST.LT.1.5D0*REPCUT(JJ)) THEN
         NNREPULSIVE=NNREPULSIVE+1
         NREPI(NNREPULSIVE)=REPI(JJ)
         NREPJ(NNREPULSIVE)=REPJ(JJ)
         NREPCUT(NNREPULSIVE)=REPCUT(JJ)
         GOTO 246
      ENDIF
   ENDDO 
 CALL FLUSH(6,ISTAT)
   DO KK=2,INTIMAGE+2 ! now check internal minima within threshold
      DMIN=1.0D10
      NI2=NOPT*(KK-2)+3*(REPI(JJ)-1)
      NI1=NOPT*(KK-1)+3*(REPI(JJ)-1)
      NJ2=NOPT*(KK-2)+3*(REPJ(JJ)-1)
      NJ1=NOPT*(KK-1)+3*(REPJ(JJ)-1)
      R1AX=XYZ(NI2+1); R1AY=XYZ(NI2+2); R1AZ=XYZ(NI2+3)
      R1BX=XYZ(NJ2+1); R1BY=XYZ(NJ2+2); R1BZ=XYZ(NJ2+3)
      R2AX=XYZ(NI1+1); R2AY=XYZ(NI1+2); R2AZ=XYZ(NI1+3)
      R2BX=XYZ(NJ1+1); R2BY=XYZ(NJ1+2); R2BZ=XYZ(NJ1+3)
      CALL INTMINONLY(R1AX,R1AY,R1AZ,R2AX,R2AY,R2AZ,R1BX,R1BY,R1BZ,R2BX,R2BY,R2BZ,DMIN,NOINT)

!     IF ((REPI(JJ).EQ.143).AND.(REPJ(JJ).EQ.191)) THEN
!        PRINT '(A,3G20.10)',' checkrep> R1AX,R1AY,R1AZ=',R1AX,R1AY,R1AZ
!        PRINT '(A,3G20.10)',' checkrep> R1BX,R1BY,R1BZ=',R1BX,R1BY,R1BZ
!        PRINT '(A,3G20.10)',' checkrep> R2AX,R2AY,R2AZ=',R2AX,R2AY,R2AZ
!        PRINT '(A,3G20.10)',' checkrep> R2BX,R2BY,R2BZ=',R2BX,R2BY,R2BZ
!        PRINT '(A,3I6,2G20.10)','JJ,REPI(JJ),REPJ(JJ),LDIST,DMIN=',JJ,REPI(JJ),REPJ(JJ),LDIST,DMIN
!     ENDIF
      IF (NOINT) CYCLE
      IF (DMIN.LT.1.5D0*REPCUT(JJ)) THEN
         NNREPULSIVE=NNREPULSIVE+1
         NREPI(NNREPULSIVE)=REPI(JJ)
         NREPJ(NNREPULSIVE)=REPJ(JJ)
         NREPCUT(NNREPULSIVE)=REPCUT(JJ)
         GOTO 246
      ENDIF
   ENDDO 
 CALL FLUSH(6,ISTAT)
246 CONTINUE
ENDDO
IF (DEBUG) PRINT '(A,2I8)',' checkrep> number of active repulsions and total=',NNREPULSIVE,NREPULSIVE

END SUBROUTINE CHECKREP

SUBROUTINE RWG(NITER,INTIMAGE,XYZ)
USE PORFUNCS
USE KEY,ONLY: FILTH,FILTHSTR,STOCKT,AMHT,SEQ,NUMGLY,STOCKAAT, RBAAT
USE COMMONS, ONLY: ZSYM, NRBSITES 
USE AMHGLOBALS, ONLY : NMRES
USE COMMONS, ONLY: NATOMS, NOPT
IMPLICIT NONE
CHARACTER(LEN=10) :: XYZFILE   = 'int.xyz   '
CHARACTER(LEN=12) :: RBXYZFILE = 'rbint.xyz   '
INTEGER,INTENT(IN) :: NITER
INTEGER :: J1,J2,GLY_COUNT,INTIMAGE
CHARACTER(LEN=80) :: FILENAME,FILENAME2,DUMMYS,DUMMYS2
DOUBLE PRECISION XYZ(NOPT*(INTIMAGE+2))

IF (FILTH.EQ.0) THEN
   FILENAME=XYZFILE
   IF (RBAAT) FILENAME2=RBXYZFILE
ELSE
   FILENAME=TRIM(XYZFILE)//'.'//TRIM(ADJUSTL(FILTHSTR))
   IF (RBAAT) FILENAME2=TRIM(RBXYZFILE)//'.'//TRIM(ADJUSTL(FILTHSTR))
ENDIF 

IF (NITER.GT.0) THEN
   IF (FILTH.EQ.0) THEN
      WRITE(DUMMYS,'(I8)') NITER
      DUMMYS2=TRIM(ADJUSTL(FILENAME))
      FILENAME='int.' // TRIM(ADJUSTL(DUMMYS)) // '.xyz' ! so that vmd recognises the file type!
      FILENAME2='rbint.' // TRIM(ADJUSTL(DUMMYS)) // '.xyz'
   ELSE 
      WRITE(DUMMYS,'(I8)') NITER
      DUMMYS2=TRIM(ADJUSTL(FILENAME))
      FILENAME='int.' // TRIM(ADJUSTL(DUMMYS)) // '.' // TRIM(ADJUSTL(FILTHSTR)) // '.xyz' 
      FILENAME2='rbint.' // TRIM(ADJUSTL(DUMMYS)) // '.' // TRIM(ADJUSTL(FILTHSTR)) // '.xyz'
   ENDIF
ENDIF
OPEN(UNIT=993,FILE=FILENAME,STATUS='replace')
IF (STOCKT .OR. STOCKAAT) THEN
   DO J2=1,INTIMAGE+2 
      WRITE(993,'(i4/)') (natoms/2)
      DO J1=1,(natoms/2) 
         WRITE(993,'(a5,1x,6f20.10)') ZSYM((j1+2)/3), &
  & XYZ((J2-1)*NOPT+3*(J1-1)+1), XYZ((J2-1)*NOPT+3*(J1-1)+2), XYZ((J2-1)*NOPT+3*(J1-1)+3), &
  &    XYZ((J2-1)*NOPT+3*((natoms/2)+J1-1)+1), XYZ((J2-1)*NOPT+3*((natoms/2)+J1-1)+2), XYZ((J2-1)*NOPT+3*((natoms/2)+J1-1)+3)
      ENDDO
   ENDDO
ELSEIF (RBAAT .AND. (.NOT. STOCKAAT)) THEN
   PRINT '(A)',' intlbfgs> ERROR *** RGW routine needs to be taught STXYZ for this potential'
   STOP
!  OPEN(UNIT=114,FILE=FILENAME2,STATUS='unknown')
!  DO J2=1,INTIMAGE+2
!     WRITE(993,'(i4/)') NATOMS/2
!     DO J1=1,(NATOMS/2) 
!        WRITE(993,'(a5,1x,3f20.10)') 'O', &
! & XYZ((J2-1)*NOPT+3*(J1-1)+1), XYZ((J2-1)*NOPT+3*(J1-1)+2), XYZ((J2-1)*NOPT+3*(J1-1)+3)
!     ENDDO
!     CALL SITEPOS(XYZ((J2-1)*NOPT+1:J2*NOPT),STXYZ)
!     WRITE(114,'(i4/)') (NATOMS/2)*NRBSITES
!     DO J1=1,(NATOMS/2)*NRBSITES
!        J3 = 3*J1
!        WRITE(114,'(a5,1x,3f20.10)') 'O', STXYZ(J3-2), STXYZ(J3-1), STXYZ(J3)
!     ENDDO
!  ENDDO
!  CLOSE(UNIT=114)
ELSEIF (AMHT) THEN
   DO J2=1,INTIMAGE+2
!  GLY set getparams.f
!               WRITE(993,'(i4)')NATOMS +NUMGLY
!  GLY printing turned off DJW
      WRITE(993,'(i4)')NATOMS
      WRITE(993,*)'Energy'
      GLY_COUNT = 0

      DO J1=1,NMRES
         IF (SEQ(J1).EQ.8) THEN
            WRITE(993,'(a5,1x,3f20.10)') 'C1   ',XYZ((J2-1)*NOPT+9*(J1-1)+1-GLY_COUNT*3),XYZ((J2-1)*NOPT+9*(J1-1)+2-GLY_COUNT*3), &
     &                                  XYZ((J2-1)*NOPT+9*(J1-1)+3-GLY_COUNT*3)
!  GLY printing turned off DJW
!           WRITE(993,'(a5,1x,3f20.10)') 'C1   ',XYZ((J2-1)*NOPT+9*(J1-1)+1-GLY_COUNT*3),XYZ((J2-1)*NOPT+9*(J1-1)+2-GLY_COUNT*3), &
!    &                                  XYZ((J2-1)*NOPT+9*(J1-1)+3-GLY_COUNT*3)
            WRITE(993,'(a5,1x,3f20.10)') 'O    ',XYZ((J2-1)*NOPT+9*(J1-1)+4-GLY_COUNT*3),XYZ((J2-1)*NOPT+9*(J1-1)+5-GLY_COUNT*3), &
     &                                  XYZ((J2-1)*NOPT+9*(J1-1)+6-GLY_COUNT*3)
            GLY_COUNT = GLY_COUNT +1
         ELSE
            WRITE(993,'(a5,1x,3f20.10)') 'C1   ',XYZ((J2-1)*NOPT+9*(J1-1)+1-GLY_COUNT*3),XYZ((J2-1)*NOPT+9*(J1-1)+2-GLY_COUNT*3), &
     &                                  XYZ((J2-1)*NOPT+9*(J1-1)+3-GLY_COUNT*3)
            WRITE(993,'(a5,1x,3f20.10)') 'C2   ',XYZ((J2-1)*NOPT+9*(J1-1)+4-GLY_COUNT*3),XYZ((J2-1)*NOPT+9*(J1-1)+5-GLY_COUNT*3), &
     &                                  XYZ((J2-1)*NOPT+9*(J1-1)+6-GLY_COUNT*3)
            WRITE(993,'(a5,1x,3f20.10)') 'O    ',XYZ((J2-1)*NOPT+9*(J1-1)+7-GLY_COUNT*3),XYZ((J2-1)*NOPT+9*(J1-1)+8-GLY_COUNT*3), &
     &                                  XYZ((J2-1)*NOPT+9*(J1-1)+9-GLY_COUNT*3)
         ENDIF
      ENDDO
   ENDDO
ELSE
   DO J2=1,INTIMAGE+2
      WRITE(993,'(i4/)') natoms
      WRITE(993,'(a5,1x,3f20.10)') (ZSYM((j1+2)/3),xyz( (j2-1)*Nopt+j1),&
    & XYZ((J2-1)*NOPT+J1+1), XYZ((J2-1)*NOPT+J1+2),J1=1,NOPT,3)
   ENDDO
ENDIF

PRINT *, 'rwg> Interpolated image coordinates were saved to xyz file "'//TRIM(FILENAME)//'"'

CLOSE(UNIT=993)
END SUBROUTINE RWG

SUBROUTINE WRITEPROFILE(NITER,EEE,INTIMAGE)
USE KEY,ONLY: FILTH,FILTHSTR
IMPLICIT NONE 
INTEGER,INTENT(IN) :: NITER, INTIMAGE
INTEGER :: I,UNIT
DOUBLE PRECISION :: EEE(INTIMAGE+2)
CHARACTER(LEN=20) :: FILENAME

UNIT=992
IF (NITER.GT.0) THEN
   WRITE(FILENAME,'(I8)') NITER
   FILENAME='int.EofS.' // TRIM(ADJUSTL(FILENAME))
ELSE   
   FILENAME='int.EofS'
ENDIF
IF (.NOT.FILTH==0) THEN
   FILENAME=TRIM(FILENAME)//'.'//TRIM(ADJUSTL(FILTHSTR))
ENDIF
OPEN(UNIT=UNIT,FILE=FILENAME,STATUS='replace')

WRITE(UNIT=UNIT,FMT='(2g24.13)') EEE(1)
DO I=2,INTIMAGE+1
   WRITE(UNIT=UNIT,FMT='(2G24.13)') EEE(I)
ENDDO
WRITE(UNIT=UNIT,FMT='(2G24.13)') EEE(INTIMAGE+2)

CLOSE(UNIT)
PRINT '(A)',' writeprofile> Interpolated energy profile was saved to file "'//trim(filename)//'"'

END SUBROUTINE WRITEPROFILE
