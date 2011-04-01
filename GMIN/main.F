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
      NPAR=1
      MYNODE=0
      MDCRD_UNIT=20000
      MDINFO_UNIT=21000
      
      LOG_FH=22979+1
      LOG_FN="GMIN.log"
      CALL OPENF(LOG_FH,">",LOG_FN)

      WRITE(LOG_FH,'(A)') "Starting serial execution" 

!op226> Allocate memory; open files; initialize different things  {{{ 

      CALL COUNTATOMS(LOG_FH)
      NATOMS = NUMBER_OF_ATOMS

      ALLOCATE(ANV(NATOMS,NATOMS,3))         

      ALLOCATE(FIN(3*NATOMS))
      ALLOCATE(XICOM(3*NATOMS),PCOM(3*NATOMS))
      ALLOCATE(SCREENC(3*NATOMS))
      ALLOCATE(IATNUM(NATOMS), VT(NATOMS), ZSYM(NATOMS))

      VT(1:NATOMS)=0.0D0 ! TO PREVENT READING FROM UNINITIALISED MEMORY

      CALL KEYWORD
      
      ALLOCATE(FF(NSAVE),QMIN(NSAVE))
      ALLOCATE(QMINP(NSAVE,3*NATOMS))

      QMINP(1:NSAVE,1:3*NATOMS)=0.0D0 ! to prevent reading from uninitialised memory
      COORDSO(1:3*NATOMS,1:NPAR)=0.0D0 ! to prevent reading from uninitialised memory
      FF(1:NSAVE)=0 ! to prevent reading from uninitialised memorY
      VATO(1:NATOMS,1:NPAR)=0.0D0 ! to prevent reading from uninitialised memory
      ALLOCATE(ESAVE(NTAB,NPAR),XINSAVE(NTAB,NPAR))
      ALLOCATE(VEC(NVEC))

      IF (PAIRDISTT) THEN
         MYPUNIT=3000+MYNODE
         CALL OPENF(PAIRDIST_FH,">>","pairdists")
         WRITE(MYPUNIT,'(A10)',ADVANCE="NO") "Quench  "
         DO J1=1,NPAIRS
            WRITE(atom1,*) PAIRDIST(J1,1)
            WRITE(atom2,*) PAIRDIST(J1,2)
            WRITE(atompair,*) TRIM(ADJUSTL(atom1))//"-"//TRIM(ADJUSTL(atom2))
            WRITE(MYPUNIT,'(A10)',ADVANCE="NO") TRIM(ADJUSTL(atompair))//"  " 
         ENDDO
         WRITE(LOG_FH,'(A)') ""
      ENDIF

      IF (TRACKDATAT) THEN
         CALL OPENF(ENERGY_FH,">>","energy")
         CALL OPENF(MARKOV_FH,">>","markov")
         CALL OPENF(BEST_FH,">>","best")
         IF (RMST) CALL OPENF(RMSD_FH,">>","rmsd")
      ENDIF
      
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

      IF (SUPERSTEP) NSUPERSTEP=0
      
      DO JP=1,NPAR
         NQ(JP)=1
      ENDDO
      DO J1=1,NSAVE
         QMIN(J1)=1.0D10
      ENDDO
! }}}

      CALL MCRUNS(SCREENC)

      IF (ALLOCATED(FIN)) DEALLOCATE(FIN)
      IF (ALLOCATED(XICOM)) DEALLOCATE(XICOM)
      IF (ALLOCATED(PCOM)) DEALLOCATE(PCOM)
      IF (ALLOCATED(ANV)) DEALLOCATE(ANV)

      CALL FLUSH(LOG_FH)

      CLOSE(LOG_FH)
      IF (PAIRDISTT) CLOSE(PAIRDIST_FH)
      IF (TRACKDATAT) THEN
         CLOSE(ENERGY_FH)
         CLOSE(MARKOV_FH)
         IF (RMST) CLOSE(RMSD_FH)
         CLOSE(BEST_FH)
      ENDIF

      STOP

      END
