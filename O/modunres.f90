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
MODULE MODUNRES
   IMPLICIT NONE
   SAVE

! JMC UNRES RELATED KEYWORDS

      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: UREFCOORD,UREFPPSANGLE ! STORE REFERENCE CARTESIANS AND INTERNALS
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: INTSTEP ! INTERNAL COORDINATE STEP VECTOR OFF TS
! THESE 2 LINES ARE FOR USE WITH CONSEC (SURPRISINGLY ENOUGH) KEYWORD TO INTERPOLATE INTERNAL COORDS OVER A SECTION OF THE 
! MOLECULE ONLY (SELECTED A/C TO RESIDUE NUMBERS AND STORED IN [START/END]RES() AS APPROPRIATE.)
      LOGICAL :: CONSECT
      INTEGER :: STARTRES(10),ENDRES(10),NUMSEC
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: MYQMINSAVE
      INTEGER COUNTER

! JMC UNRES.H
! ADDED COMMON /INTERACT/ TO KNOW ITYPE()

      INTEGER :: UNMAXRES, UNMAXRES2, MAXDIM, MAXVAR, NTYP, NTYP1
      PARAMETER (UNMAXRES=100,UNMAXRES2=2*UNMAXRES,MAXDIM=(UNMAXRES-1)*(UNMAXRES-2)/2,MAXVAR=4*UNMAXRES)
      PARAMETER (NTYP=20,NTYP1=NTYP+1)

      INTEGER :: NRES,NRES0,IALPH,IVAR,NTHETA,NPHI,NSIDE,NVARU,NFL,ICG
      DOUBLE PRECISION :: THETA,PHI,ALPH,OMEG
      DOUBLE PRECISION :: C,DC,XLOC,XROT,DC_NORM, &
     &                    DCDV,DXDV,DXDS,GRADX,GRADC,GVDWC,GELC,GRADX_SCP, &
     &                    GVDWC_SCP,GHPBX,GHPBC,GLOC,GRADCORR,GRADXORR,GVDWX

      COMMON /CHAIN/ C(3,UNMAXRES2+2),DC(3,UNMAXRES2),XLOC(3,UNMAXRES), &
     & XROT(3,UNMAXRES),DC_NORM(3,UNMAXRES2),NRES,NRES0

      DOUBLE PRECISION :: DSC,VBL,VBLINV,VBLINV2,VBL_CIS,VBL0
      COMMON /SCLOCAL/ DSC(NTYP1)
      COMMON /PEPTBOND/ VBL,VBLINV,VBLINV2,VBL_CIS,VBL0

      COMMON /VAR/ THETA(UNMAXRES),PHI(UNMAXRES),ALPH(UNMAXRES),OMEG(UNMAXRES), &
     &          IALPH(UNMAXRES,2),IVAR(4*UNMAXRES2),NTHETA,NPHI,NSIDE,NVARU

      COMMON /DERIVAT/ DCDV(6,MAXDIM),DXDV(6,MAXDIM),DXDS(6,UNMAXRES), &
     &                 GRADX(3,UNMAXRES,2),GRADC(3,UNMAXRES,2),GVDWX(3,UNMAXRES), &
     &                 GVDWC(3,UNMAXRES),GELC(3,UNMAXRES),GRADX_SCP(3,UNMAXRES), &
     &                 GVDWC_SCP(3,UNMAXRES),GHPBX(3,UNMAXRES),GHPBC(3,UNMAXRES), &
     &                 GLOC(MAXVAR,2),GRADCORR(3,UNMAXRES),GRADXORR(3,UNMAXRES),NFL,ICG

      CHARACTER(LEN=3) :: RESTYP
      CHARACTER(LEN=1) :: ONELET
      COMMON /NAMES/ RESTYP(NTYP+1),ONELET(NTYP+1)

      INTEGER :: NNT,NCT,NINT_GR,ISTART,IEND,ITYPE,MAXINT_GR,EXPON,EXPON2, &
     &           ITEL,ITYPRO,IELSTART,IELEND,NSCP_GR,ISCPSTART,ISCPEND,IATSC_S, &
     &           IATSC_E,IATEL_S,IATEL_E,IATSCP_S,IATSCP_E,ISPP,ISCP
      PARAMETER (MAXINT_GR=2)
      DOUBLE PRECISION :: AA,BB,AUGM,AAD,BAD,APPU,BPP,AEL6,AEL3

      COMMON /INTERACT/AA(NTYP,NTYP),BB(NTYP,NTYP),AUGM(NTYP,NTYP), &
     & AAD(NTYP,2),BAD(NTYP,2),APPU(2,2),BPP(2,2),AEL6(2,2),AEL3(2,2), &
     & EXPON,EXPON2,NNT,NCT,NINT_GR(UNMAXRES),ISTART(UNMAXRES,MAXINT_GR), &
     & IEND(UNMAXRES,MAXINT_GR),ITYPE(UNMAXRES),ITEL(UNMAXRES),ITYPRO, &
     & IELSTART(UNMAXRES),IELEND(UNMAXRES),NSCP_GR(UNMAXRES), &
     & ISCPSTART(UNMAXRES,MAXINT_GR),ISCPEND(UNMAXRES,MAXINT_GR), &
     & IATSC_S,IATSC_E,IATEL_S,IATEL_E,IATSCP_S,IATSCP_E,ISPP,ISCP

!     INTEGER LTM(UNMAXRES,200)
!     COMMON /RESIDINFO/ LTM

! MASS INFO FOR CALPHA, 'D' CAPPING GROUP (GLY) AND FOR THE SIDE CHAINS.
! SIDE CHAIN MASSES ARE THE SUMS OF THE ATOMIC MASSES.
      DOUBLE PRECISION :: MASSES(NTYP+2)
      COMMON /UNRESMASS/ MASSES

! VALUES SET IN INITIALIZE_P.F VIA BLOCK DATA.
! ORDER IS (FOLLOWING UNRES' RESTYP (IN INITIALIZE_P.F):
! C, CYS, MET , PHE , ILE , LEU , VAL , TRP , TYR , ALA , GLY , THR ,
! SER , GLN , ASN , GLU , ASP , HIS , ARG , LYS , PRO , D 
! DATA FROM HTTP://WWW.PHARMACY.PURDUE.EDU/~MCMP304/AATUTORIAL/
! ALL AT PHYSIOLOGICAL PH, THEREFORE ARG,ASP,GLU,LYS ARE CHARGED.
!     MASSES(1)=12.01
!     MASSES(2)=47.1
!     MASSES(3)=75.1
!     MASSES(4)=91.1
!     MASSES(5)=57.1
!     MASSES(6)=57.1
!     MASSES(7)=43.1
!     MASSES(8)=130.2
!     MASSES(9)=107.1
!     MASSES(10)=15.0
!     MASSES(11)=1.0
!     MASSES(12)=45.1
!     MASSES(13)=31.0
!     MASSES(14)=72.1
!     MASSES(15)=58.1
!     MASSES(16)=72.1
!     MASSES(17)=58.1
!     MASSES(18)=81.1
!     MASSES(19)=101.2
!     MASSES(20)=73.1
!     MASSES(21)=41.1
!     MASSES(22)=0.0

END MODULE
