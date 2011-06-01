!op226>=================================== 
!op226> GPL License info {{{
C   GMIN: A program for finding global minima
C   Copyright (C) 1999- David J. Wales
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
      
!> \name GMIN
!> \brief Main program file for GMIN
!> \param SCREENC
!op226>=================================== 
      PROGRAM GMIN
!op226>=================================== 
!op226> Declarations {{{ 
      !op226> Modules {{{ 
      USE MODMXATMS
      USE NOA
      USE COMMONS
      USE QMODULE
      USE PERMU
      USE F1COM
      USE MODAMBER
      USE MODAMBER9, only : AMBFINALIO_NODE,MDCRD_UNIT,MDINFO_UNIT,AMBPDB_UNIT
      USE MODCHARMM
      USE PORFUNCS
      !op226>  }}}
      !op226> Other {{{
      IMPLICIT NONE
      !EXTERNAL READ_CMD_ARGS
#ifdef MPI
      INCLUDE 'mpif.h'
#endif
      INTEGER J1,J2, JP, MPIERR, NDUMMY3,NPTOTAL,VERSIONTEMP
      DOUBLE PRECISION, ALLOCATABLE :: SCREENC(:)
      DOUBLE PRECISION POTEL
      DOUBLE PRECISION, ALLOCATABLE :: TMPCOORDS(:)
      INTEGER, ALLOCATABLE :: NDUMMY(:), NDUMMY2(:,:)
      LOGICAL LOPEN

      CHARACTER(LEN=130) ISTR
      CHARACTER(LEN=40) :: atom1,atom2,atompair

      COMMON /MYPOT/ POTEL
!op226>  }}}
!op226> End declarations }}}
!op226>=================================== 

      CALL CPU_TIME(TSTART)
      CALL READ_CMD_ARGS
#ifdef MPI
!op226> MPI section {{{ 
      CALL MPI_INIT(MPIERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPAR,MPIERR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYNODE,MPIERR)
      MYUNIT=22980+MYNODE
      MDCRD_UNIT=20000+MYNODE
      MDINFO_UNIT=21000+MYNODE
      AMBPDB_UNIT=18000+MYNODE
      AMBFINALIO_NODE=200+MYNODE
      WRITE (ISTR, '(I10)') MYNODE+1
      IF (TRIM(ADJUSTL(INFILE)).EQ.'') INFILE='GMIN_out'
      IF (NPAR.GT.1) THEN 
         MYFILENAME=TRIM(ADJUSTL(INFILE))//"."//TRIM(ADJUSTL(ISTR))
      ELSE
         MYFILENAME=TRIM(ADJUSTL(INFILE))
      ENDIF
      OPEN(MYUNIT,FILE=MYFILENAME, STATUS="unknown", form="formatted")
      WRITE(MYUNIT, '(A,I10,A,I10)') "Starting parallel execution: Processor", mynode+1, " of ",NPAR
!op226>}}} 
#else
!op226> Non-MPI section {{{ 
      NPAR=1
      MYNODE=0
      MYUNIT=22979+1
      MDCRD_UNIT=20000
      MDINFO_UNIT=21000
      AMBPDB_UNIT=18000
      AMBFINALIO_NODE=1
      IF (TRIM(ADJUSTL(INFILE)).EQ.'') THEN
         MYFILENAME="GMIN_out"
      ELSE
         MYFILENAME=TRIM(ADJUSTL(INFILE))
      ENDIF
      OPEN(MYUNIT,FILE=MYFILENAME, STATUS="unknown", form="formatted")
      WRITE(MYUNIT, '(A,I10,A,I10,A)') "Starting serial execution"
!op226>}}} 
#endif
!op226> Allocate memory; open files; initialize different things  {{{ 
! Add the GMIN version to the output - helps bug hunting :)
      !WRITE(MYUNIT, '(A,I5)') 'GMIN version r',VERSIONTEMP
!     CALL DISPLAY_VERSION(MYUNIT)
      CALL COUNTATOMS(MYUNIT)
      IF (CHRMMT) WRITE(MYUNIT,'(A,I8)') 'main> MAXAIM parameter for CHARMM MXATMS=',MXATMS
      CALL MODCOMMONINIT
 
      ALLOCATE(FIN(3*NATOMS))
      ALLOCATE(XICOM(3*NATOMS),PCOM(3*NATOMS))
      ALLOCATE(SCREENC(3*NATOMS))
      ALLOCATE(IATNUM(NATOMS), VT(NATOMS), ZSYM(NATOMS))
      VT(1:NATOMS)=0.0D0 ! TO PREVENT READING FROM UNINITIALISED MEMORY
      IF (Q4T) CALL SHINIT

      INQUIRE(UNIT=1,OPENED=LOPEN)
      IF (LOPEN) THEN
         WRITE(*,'(A,I2,A)') 'main> A ERROR *** Unit ', 1, ' is not free '
         STOP
      ENDIF

      CALL KEYWORD

      INQUIRE(UNIT=1,OPENED=LOPEN)
      IF (LOPEN) THEN
         WRITE(*,'(A,I2,A)') 'main> B ERROR *** Unit ', 1, ' is not free '
         STOP
      ENDIF

      
      IF (CHRIGIDTRANST.OR.CHRIGIDROTT) CALL CHSETSEG

!op226> RMST {{{ 
      IF (RMST) THEN
         ALLOCATE(RMSBEST(RMSSAVE,2),RMSCOOR(RMSSAVE,3*NATOMS))
         RMSBEST(1:RMSSAVE,1)=RMSLIMIT+RMSTOL
         RMSBEST(1:RMSSAVE,2)=0.D0
         RMSCOOR(1:RMSSAVE,1:3*NATOMS)=0.D0
         ALLOCATE(COORCOMP(1:3*NATOMS))
!
!        csw34> Need to add the read in for other file formats
!        (crd,pdb). Still not sure how GMIN reads in the input.crd file
!        for CHARMM. It appears to be reading into unit 7 in io1.f
!        though - but only from xyz format! Is there a hidden
!        conversion?
!         
         OPEN(UNIT=1,FILE='compare',STATUS='OLD')
         READ(1,*) (COORCOMP(J1),J1=1,3*NATOMS)
         CLOSE(1)
      ENDIF
!op226>}}} 

      ALLOCATE(FF(NSAVE),QMIN(NSAVE))
      ALLOCATE(QMINP(NSAVE,3*NATOMS))

!        csw34> ALLOCATE the interaction energy tracking arrays if A9INTE in data
!op226> IF (A9INTET.AND.AMBERT) THEN  {{{
      IF (A9INTET.AND.AMBERT) THEN 
         ALLOCATE(INTEFF(NSAVEINTE),INTEQMIN(NSAVEINTE))
         ALLOCATE(INTEQMINP(NSAVEINTE,3*NATOMS))
         INTEQMIN(:)=1.0D10 
         INTEQMINP(1:NSAVEINTE,1:3*NATOMS)=0.0D0 ! to prevent reading from uninitialised memory
         INTEFF(1:NSAVEINTE)=0 ! to prevent reading from uninitialised memory
      ENDIF  
!op226>}}} 

!op226> GAUSST {{{ 
      IF (GAUSST) THEN
         ALLOCATE(GAUSSKK(3*NATOMS,GMODES),GAUSSEE(GMODES),GKSMALL(3*NATOMS))
         CALL KEGEN ! INITIAL SETUP 
         DO J1=1,GMODES
             PRINT *,J1,GAUSSEE(J1)
         ENDDO
      ENDIF
!op226>}}} 

      QMINP(1:NSAVE,1:3*NATOMS)=0.0D0 ! to prevent reading from uninitialised memory
      COORDSO(1:3*NATOMS,1:NPAR)=0.0D0 ! to prevent reading from uninitialised memory
      FF(1:NSAVE)=0 ! to prevent reading from uninitialised memorY
      VATO(1:NATOMS,1:NPAR)=0.0D0 ! to prevent reading from uninitialised memory
      ALLOCATE(ESAVE(NTAB,NPAR),XINSAVE(NTAB,NPAR))
      ALLOCATE(VEC(NVEC))

!     IF (SYMMETRIZE.AND.(.NOT.CENT)) THEN
!        PRINT '(A)','Probable input error - SYMMETRIZE true but CENT false'
!        STOP
!     ENDIF

!op226> DUMPT {{{ 
      IF (DUMPT) THEN
         IF (CHRMMT) THEN
            OPEN(UNIT=719,FILE='dump.crd',STATUS='UNKNOWN')
            OPEN(UNIT=720,FILE='dump.pdb',STATUS='UNKNOWN')
         ENDIF
         DUMPXYZUNIT=40
         DUMPVUNIT=39
         DO J1=1,NPAR
            WRITE (ISTR,'(A5,I1,A4)') 'dump.',J1,'.xyz'
            J2=DUMPXYZUNIT+J1
            OPEN(UNIT=J2,FILE=TRIM(ADJUSTL(ISTR)),STATUS='UNKNOWN')
            WRITE (ISTR,'(A5,I1,A2)') 'dump.',J1,'.V'
            J2=DUMPVUNIT-J1
            OPEN(UNIT=J2,FILE=TRIM(ADJUSTL(ISTR)),STATUS='UNKNOWN')
         ENDDO
      ENDIF
!op226>}}} 

      IF (PAIRDISTT) THEN
         MYPUNIT=3000+MYNODE
         IF (NPAR.GT.1) THEN
            OPEN(MYPUNIT,FILE="pairdists."//trim(adjustl(istr)),STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
         ELSE
            OPEN(MYPUNIT,FILE="pairdists",STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
         ENDIF
         WRITE(MYPUNIT,'(A10)',ADVANCE="NO") "Quench  "
         DO J1=1,NPAIRS
            WRITE(atom1,*) PAIRDIST(J1,1)
            WRITE(atom2,*) PAIRDIST(J1,2)
            WRITE(atompair,*) TRIM(ADJUSTL(atom1))//"-"//TRIM(ADJUSTL(atom2))
            WRITE(MYPUNIT,'(A10)',ADVANCE="NO") TRIM(ADJUSTL(atompair))//"  " 
         ENDDO
         WRITE(MYUNIT,'(A)') ""
      ENDIF

!op226> TRACKDATAT {{{ 
!
!     csw34> TRACKDATA keyword prints the energy and markov energy 
!     to files for viewing during a run. If RMS is also specified it
!     prints the rmsd from the comparison structure into a file.
!
      IF (TRACKDATAT) THEN
         MYEUNIT=4000+MYNODE
         MYMUNIT=6000+MYNODE
         MYRUNIT=8000+MYNODE
         MYBUNIT=10000+MYNODE
         IF (NPAR.GT.1) THEN
            OPEN(MYEUNIT,FILE="energy."//trim(adjustl(istr)),STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
            OPEN(MYMUNIT,FILE="markov."//trim(adjustl(istr)),STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
            OPEN(MYBUNIT,FILE="best."//trim(adjustl(istr)),STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
         ELSE
            OPEN(MYEUNIT,FILE='energy',STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
            OPEN(MYMUNIT,FILE='markov',STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
            IF (RMST) OPEN(MYRUNIT,FILE='rmsd',STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
            OPEN(MYBUNIT,FILE='best',STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
            IF (A9INTET) THEN
               OPEN(UNIT=3998,FILE='intE.dat',STATUS='UNKNOWN',FORM='FORMATTED')
               OPEN(UNIT=3999,FILE='bestintE.dat',STATUS='UNKNOWN',FORM='FORMATTED')
            ENDIF            
         ENDIF
      ENDIF
!op226>}}} 
      CALL FLUSH(6)
      CALL IO1
!op226> CSMT {{{ 
!
! If this is a CSM optimisation we now have to multiply the number of atoms by the number of
! group operations and replicate some coordinates and allowed permutations.
!
      IF (CSMT) THEN
         CALL CSMINIT
         IF (SYMMETRIZECSM) THEN
            IF (CSMMAXIT.EQ.0) CSMMAXIT=MAXIT
         ELSE
            CSMMAXIT=MAXIT
         ENDIF
      ENDIF
!op226>}}} 

      CALL FLUSH(6)

!op226> PROJIT {{{ 
      IF (PROJIT) THEN
         ALLOCATE(TMPCOORDS(3*NATOMS))
         TMPCOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS,1)
         CALL PROJIINIT(TMPCOORDS,NATOMS)
         DEALLOCATE(TMPCOORDS)
      ENDIF
!op226>}}} 

!op226> VGW SEEDT {{{ 
      IF(VGW) THEN
        CALL INITIALIZE_VGWSP(NATOMS,LJSIGMA,LJEPSILON,TAUMAX,CPFACTORSG,CPS,VGWTOL)
        CALL INITIALIZE_VGW(NATOMS,LJSIGMA,LJEPSILON,TAUMAXFULL,CPFACTORFG,CPF)
      ENDIF

      IF (SEEDT) THEN
         CALL GSEED
      ELSE
         IF ((.NOT.FIELDT).AND.CENT) THEN
            DO J1=1,NPAR
               IF (.NOT.SEEDT) CALL CENTRE2(COORDS(1:3*NATOMS,J1))
            ENDDO
         ELSEIF ((.NOT.FIELDT).AND.FIXCOM) THEN
            DO J1=1,NPAR
               IF (.NOT.SEEDT) CALL CENTRECOM(COORDS(1:3*NATOMS,J1))
            ENDDO
         ENDIF
      ENDIF
!op226>}}} 

      IF (SUPERSTEP) NSUPERSTEP=0
      
      DO JP=1,NPAR
         NQ(JP)=1
      ENDDO
      DO J1=1,NSAVE
         QMIN(J1)=1.0D10
      ENDDO
!op226> End initializations and allocations }}} 

!op226> Main program run
      IF ((NRUNS.GT.0).OR.PTMC.OR.BSPT) CALL MCRUNS(SCREENC)

C     CALL SYSTEM('rm ssdump ssave >& /dev/null')

!op226> Deallocate memory; close files {{{ 
      IF (ALLOCATED(FIN)) DEALLOCATE(FIN)
      IF (ALLOCATED(XICOM)) DEALLOCATE(XICOM)
      IF (ALLOCATED(PCOM)) DEALLOCATE(PCOM)
      IF (ALLOCATED(GAUSSKK)) DEALLOCATE(GAUSSKK,GAUSSEE)
C     deallocate dynamic memory for AMBER
      IF (AMBERT) CALL AMBER_DEALLOCATE_STACKS
      CALL MODCOMMONDEINIT
      CALL FLUSH(MYUNIT)

      IF(VGW) THEN
        CALL CLEANUP_VGWSP
        CALL CLEANUP_VGW
      ENDIF

      CLOSE(MYUNIT)
! csw34> close pairdists.* files
      IF (PAIRDISTT) CLOSE(MYPUNIT)
!op226> closing files for: TRACKDATAT RMST A9INTET {{{ 
      IF (TRACKDATAT) THEN
         CLOSE(MYEUNIT)
         CLOSE(MYMUNIT)
         IF (RMST) CLOSE(MYRUNIT)
         CLOSE(MYBUNIT)
         IF (A9INTET) THEN
            CLOSE(3998)
            CLOSE(3999)
         ENDIF
      ENDIF
!op226>}}} 
#ifdef MPI
      CALL MPI_FINALIZE(MPIERR)
#endif

!op226>}}} 
      STOP

      END
