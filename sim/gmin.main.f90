      
!> @name GMIN
!
!> @brief Main program file for GMIN
!

      PROGRAM GMIN
! declarations {{{
      USE GMIN_COMMONS
      USE GMIN_PORFUNCS
      USE GMIN_FUNC
      
      IMPLICIT NONE
      
      INTEGER :: NUMBER_OF_ATOMS

      INTEGER,ALLOCATABLE :: FF(:),INTEFF(:) ! NSAVE

      DOUBLE PRECISION, ALLOCATABLE :: QMIN(:), INTEQMIN(:) ! NSAVE
      DOUBLE PRECISION, ALLOCATABLE :: QMINP(:,:), INTEQMINP(:,:)
    
      DOUBLE PRECISION, ALLOCATABLE :: FIN(:)

      DOUBLE PRECISION, ALLOCATABLE :: SCREENC(:,:)

      DOUBLE PRECISION POTEL
      DOUBLE PRECISION, ALLOCATABLE :: TMPCOORDS(:)
      INTEGER, ALLOCATABLE :: NDUMMY(:), NDUMMY2(:,:)

      ! USED IN WRITING TO PAIRDIST      
      CHARACTER(LEN=40) :: ATOM1,ATOM2,ATOMPAIR

      COMMON /MYPOT/ POTEL
! }}}

! init {{{
      CALL CPU_TIME(TSTART)
      CALL KEYWORD
      CALL READ_CMD_ARGS

      CALL OPENF(LFH,">",LFN)

      WRITE(LFH,'(A)') "Starting serial execution" 

      ALLOCATE(ANV(NATOMS,NATOMS,3))         

      ALLOCATE(FIN(NATOMS,3))
      ALLOCATE(XICOM(NATOMS,3),PCOM(NATOMS,3))
      ALLOCATE(SCREENC(NATOMS,3))
      ALLOCATE(IATNUM(NATOMS), VT(NATOMS), ZSYM(NATOMS))

      VT(1:NATOMS)=0.0D0 ! TO PREVENT READING FROM UNINITIALISED MEMORY

      
      ALLOCATE(FF(NSAVE),QMIN(NSAVE))
      ALLOCATE(QMINP(NSAVE,3*NATOMS))

      QMINP(1:NSAVE,1:3*NATOMS)=0.0D0 ! to prevent reading from uninitialised memory
      COORDSO=0.0D0 ! to prevent reading from uninitialised memory
      FF=0 ! to prevent reading from uninitialised memorY
      VATO=0.0D0 ! to prevent reading from uninitialised memory
      ALLOCATE(ESAVE(NTAB,NPAR),XINSAVE(NTAB,NPAR))
      ALLOCATE(VEC(NVEC))

      include main.pairdist.inc.f90
      
      IF (TRACKDATAT) THEN
         CALL OPENF(ENERGY_FH,">>","energy")
         CALL OPENF(MARKOV_FH,">>","markov")
         CALL OPENF(BEST_FH,">>","best")
         IF (RMST) CALL OPENF(RMSD_FH,">>","rmsd")
      ENDIF
      
      CALL FLUSH(6)
      CALL IO1
      CALL FLUSH(6)

      CALL CENTRE2(COORDS)

      NQ(1)=1

      DO J1=1,NSAVE
         QMIN(J1)=1.0D10
      ENDDO
! }}}
      CALL MCRUNS(SCREENC)
! end {{{
      IF (ALLOCATED(FIN)) DEALLOCATE(FIN)
      IF (ALLOCATED(XICOM)) DEALLOCATE(XICOM)
      IF (ALLOCATED(PCOM)) DEALLOCATE(PCOM)
      IF (ALLOCATED(ANV)) DEALLOCATE(ANV)

      CALL FLUSH(LFH)
      CLOSE(LFH)

      IF (PAIRDISTT) CLOSE(PAIRDIST_FH)
      IF (TRACKDATAT) THEN
         CLOSE(ENERGY_FH)
         CLOSE(MARKOV_FH)
         IF (RMST) CLOSE(RMSD_FH)
         CLOSE(BEST_FH)
      ENDIF
! }}}

      STOP

      END
