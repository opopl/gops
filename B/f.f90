      MODULE F

      USE V 
      USE COMMONS

      IMPLICIT NONE 

      ! interfaces  {{{
      INTERFACE

!      ! ep46 {{{
        !SUBROUTINE EP46(FH,DEB,QO, N, GRAD, E, GTEST)
            !IMPLICIT NONE
            !INTEGER,INTENT(IN) :: N
            !LOGICAL,INTENT(IN) :: DEB
            !INTEGER,INTENT(IN) :: FH
            !DOUBLE PRECISION,DIMENSION(3*N),INTENT(IN) :: QO
            !DOUBLE PRECISION,DIMENSION(3*N),INTENT(OUT) :: GRAD
            !DOUBLE PRECISION,DIMENSION(10) :: E
            !LOGICAL,INTENT(IN) :: GTEST
        !ENDSUBROUTINE EP46
        !! }}}
        !! eg46 {{{
        !SUBROUTINE EG46(FH,DEB,QO, N, GRAD, E, GTEST)
			!IMPLICIT NONE
			!INTEGER,INTENT(IN) :: N
			!LOGICAL,INTENT(IN) :: DEB
			!INTEGER,INTENT(IN) :: FH
			!DOUBLE PRECISION,DIMENSION(3*N),INTENT(IN) :: QO
			!DOUBLE PRECISION,DIMENSION(3*N),INTENT(OUT) :: GRAD
            !DOUBLE PRECISION,DIMENSION(10) :: E
			!LOGICAL,INTENT(IN) :: GTEST
        !ENDSUBROUTINE EG46
        !! }}}
        !! calc_int_coords {{{

        !SUBROUTINE CALC_INT_COORDS(QO,N,A_PARAM,B_PARAM,C_PARAM,D_PARAM,&
                !& X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD, &
     !&                             BOND_ANGLE,TOR_ANGLE,RADII,NTYPE,PTYPE)
            !IMPLICIT NONE
    
            !! sub
            !DOUBLE PRECISION,DIMENSION(3*N),INTENT(IN) :: QO
            !INTEGER,INTENT(IN) :: N 
            !DOUBLE PRECISION,DIMENSION(N,N),INTENT(IN) :: A_PARAM,B_PARAM
            !DOUBLE PRECISION,DIMENSION(N),INTENT(IN) :: C_PARAM,D_PARAM
            !DOUBLE PRECISION,DIMENSION(N),INTENT(OUT) :: X,Y,Z
            !DOUBLE PRECISION,DIMENSION(N,N),INTENT(OUT) :: XR,YR,ZR
            !DOUBLE PRECISION,DIMENSION(N,3),INTENT(OUT) :: DOT_PROD
            !DOUBLE PRECISION,DIMENSION(N),INTENT(OUT) :: X_PROD, BOND_ANGLE, TOR_ANGLE
            !DOUBLE PRECISION,DIMENSION(N,N),INTENT(OUT) :: RADII 
            !INTEGER,DIMENSION(N),INTENT(IN) :: NTYPE
            !CHARACTER(LEN=*),INTENT(IN) :: PTYPE

        !ENDSUBROUTINE CALC_INT_COORDS       

        !!  }}}
        !! calc_energy {{{
        !SUBROUTINE CALC_ENERGY(FH,DEB,QO,E,N,&
                !& A_PARAM,B_PARAM,C_PARAM,D_PARAM,&
                !& X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD, &
                !& BOND_ANGLE,TOR_ANGLE,RADII,NTYPE,PTYPE)

            !IMPLICIT NONE
            !LOGICAL DEB
            !INTEGER,INTENT(IN) :: FH
            !DOUBLE PRECISION,DIMENSION(10) :: E
            !DOUBLE PRECISION,DIMENSION(3*N),INTENT(IN) :: QO
            !INTEGER,INTENT(IN) :: N
            !DOUBLE PRECISION,DIMENSION(N,N),INTENT(IN) :: A_PARAM,B_PARAM
            !DOUBLE PRECISION,DIMENSION(N),INTENT(IN) :: C_PARAM,D_PARAM
            !DOUBLE PRECISION,DIMENSION(N),INTENT(IN) :: X,Y,Z
            !DOUBLE PRECISION,DIMENSION(N,N),INTENT(IN) :: XR,YR,ZR
            !DOUBLE PRECISION,DIMENSION(N,3),INTENT(IN) :: DOT_PROD
            !DOUBLE PRECISION,DIMENSION(N),INTENT(IN) :: X_PROD, BOND_ANGLE, TOR_ANGLE
            !DOUBLE PRECISION,DIMENSION(N,N),INTENT(IN) :: RADII 
            !INTEGER,DIMENSION(N),INTENT(IN) :: NTYPE
            !CHARACTER(LEN=*),INTENT(IN) :: PTYPE

        !ENDSUBROUTINE CALC_ENERGY
        !! }}}
        !! calc_gradient {{{
         !SUBROUTINE CALC_GRADIENT(FH,DEB,QO,FQ,N,&
               !&  A_PARAM,B_PARAM,C_PARAM,D_PARAM,X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD, &
               !&  BOND_ANGLE,TOR_ANGLE,RADII,NTYPE,PTYPE)
           !IMPLICIT NONE
    
            !LOGICAL DEB
            !DOUBLE PRECISION,DIMENSION(3*N),INTENT(IN) :: QO
            !INTEGER,INTENT(IN) :: N,FH
            !DOUBLE PRECISION,DIMENSION(N,N),INTENT(IN) :: A_PARAM,B_PARAM
            !DOUBLE PRECISION,DIMENSION(N),INTENT(IN) :: C_PARAM,D_PARAM
            !DOUBLE PRECISION,DIMENSION(N),INTENT(IN) :: X,Y,Z
            !DOUBLE PRECISION,DIMENSION(N,N),INTENT(IN) :: XR,YR,ZR
            !DOUBLE PRECISION,DIMENSION(N,3),INTENT(IN) :: DOT_PROD
            !DOUBLE PRECISION,DIMENSION(N),INTENT(IN) :: X_PROD, BOND_ANGLE, TOR_ANGLE
            !DOUBLE PRECISION,DIMENSION(N,N),INTENT(IN) :: RADII 
            !INTEGER,DIMENSION(N),INTENT(IN) :: NTYPE
            !CHARACTER(LEN=*),INTENT(IN) :: PTYPE
            !! FQ - gradient
            !DOUBLE PRECISION,DIMENSION(3*N),INTENT(OUT) :: FQ
        !ENDSUBROUTINE CALC_GRADIENT
            !! }}}
            !! calc_dyn {{{
        !SUBROUTINE CALC_DYN(QO,N,A_PARAM,B_PARAM,C_PARAM,D_PARAM,X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD, &
     !&                      BOND_ANGLE,TOR_ANGLE,RADII,NTYPE,PTYPE)
            !IMPLICIT NONE
    
            !!LOGICAL DEB
            !!INTEGER,INTENT(IN) :: FH
            !DOUBLE PRECISION,DIMENSION(3*N),INTENT(INOUT) :: QO
            !INTEGER,INTENT(IN) :: N
            !DOUBLE PRECISION,DIMENSION(N,N),INTENT(IN) :: A_PARAM,B_PARAM
            !DOUBLE PRECISION,DIMENSION(N),INTENT(IN) :: C_PARAM,D_PARAM
            !DOUBLE PRECISION,DIMENSION(N),INTENT(IN) :: X,Y,Z
            !DOUBLE PRECISION,DIMENSION(N,N),INTENT(IN) :: XR,YR,ZR
            !DOUBLE PRECISION,DIMENSION(N,3),INTENT(IN) :: DOT_PROD
            !DOUBLE PRECISION,DIMENSION(N),INTENT(IN) :: X_PROD, BOND_ANGLE, TOR_ANGLE
            !DOUBLE PRECISION,DIMENSION(N,N),INTENT(IN) :: RADII 
            !INTEGER,DIMENSION(N),INTENT(IN) :: NTYPE
            !CHARACTER(LEN=*),INTENT(IN) :: PTYPE

        !ENDSUBROUTINE CALC_DYN
            !! }}}
!            ! param_array {{{
        !SUBROUTINE PARAM_ARRAY(N,A_PARAM,B_PARAM,C_PARAM,D_PARAM,NTYPE)
			!INTEGER,INTENT(IN) :: N
			!INTEGER,DIMENSION(N),INTENT(OUT) :: NTYPE
			!DOUBLE PRECISION,DIMENSION(N,N),INTENT(OUT) :: A_PARAM,B_PARAM
			!DOUBLE PRECISION,DIMENSION(N),INTENT(OUT) :: C_PARAM,D_PARAM
        !ENDSUBROUTINE PARAM_ARRAY
        !! }}}
            !! gparam_array {{{
        !SUBROUTINE GPARAM_ARRAY(N,A_PARAM,B_PARAM,C_PARAM,D_PARAM,NTYPE)
            !implicit none
			!INTEGER,INTENT(IN) :: N
			!INTEGER,DIMENSION(N),INTENT(OUT) :: NTYPE
			!DOUBLE PRECISION,DIMENSION(N,N),INTENT(OUT) :: A_PARAM,B_PARAM
			!DOUBLE PRECISION,DIMENSION(N),INTENT(OUT) :: C_PARAM,D_PARAM
        !ENDSUBROUTINE GPARAM_ARRAY
        !! }}}
        !dumpstate {{{
        SUBROUTINE DUMPSTATE(NDONE,EBEST,BESTCOORDS,JBEST,JP)
            USE COMMONS
            USE V
            IMPLICIT NONE
            INTEGER,INTENT(IN) :: NDONE, JBEST(NPAR), JP
            DOUBLE PRECISION,INTENT(IN) ::  EBEST(NPAR), BESTCOORDS(3*NATOMS,NPAR)
        ENDSUBROUTINE
        ! }}}
        ! transition accrej mc {{{
        SUBROUTINE TRANSITION(ENEW,EOLD,ATEST,NP,RANDOM,MCTEMP)
		      DOUBLE PRECISION,INTENT(IN) :: ENEW, EOLD, MCTEMP
		      DOUBLE PRECISION,INTENT(OUT) :: RANDOM
		      LOGICAL,INTENT(OUT) :: ATEST
		      INTEGER,INTENT(IN) :: NP
        END SUBROUTINE TRANSITION

        SUBROUTINE ACCREJ(NSUCCESS,NFAIL,JP,NSUCCESST,NFAILT)
              INTEGER,DIMENSION(:),INTENT(INOUT) ::  NSUCCESS, NFAIL, NFAILT, NSUCCESST
              INTEGER,INTENT(IN) :: JP
        END SUBROUTINE ACCREJ

        SUBROUTINE MC(NSTEPS,SCALEFAC,SCREENC)
	      INTEGER ::    NSTEPS
	      DOUBLE PRECISION ::   SCALEFAC
	      DOUBLE PRECISION, DIMENSION(:) ::   SCREENC
        END SUBROUTINE MC 

        SUBROUTINE QUENCH(QTEST,NP,ITER,TIME,BRUN,QDONE,P)

	      DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: P   
	      INTEGER NP,ITER,BRUN,QDONE
	      LOGICAL QTEST
	      DOUBLE PRECISION ::   TIME

        END SUBROUTINE QUENCH

        ! }}}

        SUBROUTINE FINALQ
        ENDSUBROUTINE FINALQ

        SUBROUTINE MYRESET(JP,NATOMS,NPAR,NSEED)
            INTEGER JP,NATOMS,NPAR,NSEED
        END SUBROUTINE MYRESET

        !SUBROUTINE GSAVEIT(EREAL,P,NP)
            !INTEGER,intent(in) :: NP
            !DOUBLE PRECISION,intent(in) :: EREAL
            !DOUBLE PRECISION,intent(in), DIMENSION(:) ::   P
        !END SUBROUTINE GSAVEIT

        SUBROUTINE MYLBFGS(N,M,XCOORDS,DIAGCO,EPS,MFLAG,ENERGY,ITMAX,ITDONE,RESET,NP)

	      INTEGER :: N,M,ITMAX,ITDONE,NP
	      DOUBLE PRECISION,DIMENSION(:) ::   XCOORDS
	      LOGICAL DIAGCO,MFLAG,RESET
	      DOUBLE PRECISION ::   EPS, ENERGY

        END SUBROUTINE MYLBFGS

        
      SUBROUTINE POTENTIAL(X,GRAD,EREAL,GRADT,SECT)
	      DOUBLE PRECISION, DIMENSION(:),INTENT(IN) :: X
	      DOUBLE PRECISION, DIMENSION(:),INTENT(OUT) :: GRAD
	      DOUBLE PRECISION, INTENT(OUT) :: EREAL
	      LOGICAL, INTENT(IN) :: GRADT, SECT       
      END SUBROUTINE POTENTIAL

      END INTERFACE
      ! }}}

      CONTAINS

      ! string subs {{{
      subroutine wd(f,s,num)

      integer,intent(in) :: f,num
      character(len=*),intent(in) :: s
      character(len=100) sl
      integer i

      sl=s
      do i=1,num
        sl=sl // s
      enddo
      write(f,'(a)') sl

      endsubroutine wd

      subroutine ed(f)
      integer f
      write(f,'(a)') "==========================================="
      endsubroutine ed
      ! }}}
      ! check_file {{{
	   FUNCTION CHECK_FILE(FILE_UNIT, FOR_READ, FOR_WRITE)
	   IMPLICIT NONE
	   ! Function
	   LOGICAL             :: CHECK_FILE
	   ! Arguments
	   INTEGER, INTENT(IN) :: FILE_UNIT
	   LOGICAL, INTENT(IN) :: FOR_READ, FOR_WRITE
	   ! Local variables
	   CHARACTER (LEN=7)   :: FILE_READABLE, FILE_WRITABLE
	   LOGICAL             :: FILE_EXISTS, FILE_OPEN
	   
	   ! Initialise to .FALSE. (i.e. file isn't good for use)
	   CHECK_FILE = .FALSE.
	   
	   INQUIRE(FILE_UNIT, EXIST=FILE_EXISTS, OPENED=FILE_OPEN, READ=FILE_READABLE, WRITE=FILE_WRITABLE)
	   ! Did the file open correctly and connect to the unit?
	   IF (.NOT. (FILE_EXISTS .AND. FILE_OPEN)) THEN
	      RETURN
	   END IF
	   ! Can we read, if read access has been requested?
	   IF (FOR_READ) THEN
	      IF (FILE_READABLE .NE. 'YES') THEN
	         RETURN
	      END IF
	   END IF
	   ! Can we write, if write access has been requested?
	   IF (FOR_WRITE) THEN
	      IF (FILE_WRITABLE .NE. 'YES') THEN
	         RETURN
	      END IF
	   END IF
	   ! If all the other checks succeed, set CHECK_FILE to .TRUE., the file is ok
	   CHECK_FILE = .TRUE.
	   
	   RETURN
	   END FUNCTION CHECK_FILE
       ! }}}
      ! write_coords  {{{
      SUBROUTINE WRITE_COORDS(FILE_UNIT, FORMAT_SPEC, RUN_NUMBER)
         ! Does a sanity check and then writes COORDS to the specified unit with a 
         ! specified format.  Optional check for run number if we're doing parallel stuff.
         ! Commons
         USE COMMONS, ONLY : COORDS, NPAR
         IMPLICIT NONE
         ! Arguments
         INTEGER, INTENT(IN)           :: FILE_UNIT
         CHARACTER (LEN=*), INTENT(IN) :: FORMAT_SPEC
         ! Optional arguments
         INTEGER, INTENT(IN), OPTIONAL :: RUN_NUMBER
         ! Local variables
         INTEGER                       :: COUNTER
         
         ! Some quick sanity checks to make sure that input makes sense
         IF (PRESENT(RUN_NUMBER)) THEN
            IF (RUN_NUMBER .GT. NPAR) THEN
               STOP 'The run number is larger than the number of parallel runs.  Cannot write the output file.'
            END IF
         END IF
         
         ! Sanity checks on the file that we're writing to
         IF (.NOT. CHECK_FILE(FILE_UNIT, FOR_READ=.FALSE., FOR_WRITE=.TRUE.)) THEN
            STOP 'File did not open correctly.'
         END IF
         
         ! Write coords
         IF (PRESENT(RUN_NUMBER)) THEN
            WRITE(FILE_UNIT, FORMAT_SPEC) COORDS(:, RUN_NUMBER)
         ELSE
            WRITE(FILE_UNIT, FORMAT_SPEC) COORDS(:, 1)
         END IF
         
      END SUBROUTINE WRITE_COORDS
      ! }}}
      ! write_markov_coords {{{
      
      SUBROUTINE WRITE_MARKOV_COORDS(FILE_UNIT, FORMAT_SPEC, RUN_NUMBER)
         ! Does a sanity check and then writes COORDSO (the Markov coords) to the specified 
         ! unit with a specified format.  Optional check for run number if we're doing 
         ! parallel stuff.
         ! Commons
         USE COMMONS, ONLY : COORDSO, NPAR
         IMPLICIT NONE
         ! Arguments
         INTEGER, INTENT(IN)           :: FILE_UNIT
         CHARACTER (LEN=*), INTENT(IN) :: FORMAT_SPEC
         ! Optional arguments
         INTEGER, INTENT(IN), OPTIONAL :: RUN_NUMBER
         ! Local variables
         INTEGER                       :: COUNTER
      
         ! Some quick sanity checks to make sure that input makes sense
         IF (PRESENT(RUN_NUMBER)) THEN
            IF (RUN_NUMBER .GT. NPAR) THEN
               STOP 'The run number is larger than the number of parallel runs.  Cannot write the output file.'
            END IF
         END IF
         
         ! Sanity checks on the file that we're writing to
         IF (.NOT. CHECK_FILE(FILE_UNIT, FOR_READ=.FALSE., FOR_WRITE=.TRUE.)) THEN
            STOP 'File did not open correctly.'
         END IF
         
         ! Write coordso
         IF (PRESENT(RUN_NUMBER)) THEN
            WRITE(FILE_UNIT, FORMAT_SPEC) COORDSO(:, RUN_NUMBER)
         ELSE
            WRITE(FILE_UNIT, FORMAT_SPEC) COORDSO(:, 1)
         END IF
         
      END SUBROUTINE WRITE_MARKOV_COORDS
      ! }}}
      ! gsaveit {{{
      SUBROUTINE GSAVEIT(EREAL,P,NP)

      USE COMMONS
      USE V

      IMPLICIT NONE
      ! sub 
      INTEGER NP
      DOUBLE PRECISION :: EREAL
      DOUBLE PRECISION, DIMENSION(:) ::   P
      ! loc

      INTEGER J1, J2, J3,  NQTOT, NPCALL, CSMIT
      DOUBLE PRECISION AVVAL, CSMRMS

      COMMON /TOT/ NQTOT
      ! number of function calls
      COMMON /PCALL/ NPCALL
!
!  Save the lowest NSAVE distinguishable configurations.
!
!     WRITE(*,'(A,12E15.7)') 'EREAL,ECONV,QMIN',EREAL,ECONV,(QMIN(J1),J1=1,NSAVE)
      DO J1=1,NSAVE
         IF (DABS(EREAL-QMIN(J1)).LT.ECONV) THEN
!
!  These are probably the same - but just to make sure we save the 
!  lowest.
!
            IF (EREAL.LT.QMIN(J1)) THEN
               QMIN(J1)=EREAL
               DO J2=1,3*NATOMS
                  QMINP(J1,J2)=P(J2)
               ENDDO
            ENDIF
            GOTO 10
         ENDIF
         IF (EREAL.LT.QMIN(J1)) THEN

            J2=NSAVE
20          CONTINUE
      
            IF (NSAVE.GT.1) THEN
               QMIN(J2)=QMIN(J2-1)
               FF(J2)=FF(J2-1)
               DO J3=1,3*NATOMS
                  QMINP(J2,J3)=QMINP(J2-1,J3)
               ENDDO

               J2=J2-1
               IF (J2.GE.J1+1) GOTO 20
            ENDIF

            QMIN(J1)=EREAL
            FF(J1)=NQ(NP)
            DO J2=1,3*NATOMS
               QMINP(J1,J2)=P(J2)
            ENDDO

            GOTO 10
         ENDIF
      ENDDO

10    CONTINUE

      DO J1=1,NTARGETS
         IF (EREAL-TARGETS(J1).LT.ECONV) THEN
            WRITE(LFH,'(2(A,I15),A)') 'saveit> Target hit after ',NQTOT,' quenches ',NPCALL,' function calls.'
            WRITE(LFH,'(2(A,F20.10))') 'saveit> Energy=',EREAL,' target=',TARGETS(J1)
            HIT=.TRUE.
         ENDIF
      ENDDO

      RETURN
      ENDSUBROUTINE GSAVEIT
      ! }}}
      ! trans gseed reseed pairdistance  {{{

      DOUBLE PRECISION FUNCTION TRANS(X,XMIN,GAMMA)
      IMPLICIT NONE
      DOUBLE PRECISION X, XMIN, GAMMA

      TRANS=1.0D0 - EXP(-GAMMA*(X-XMIN))

      RETURN
      END FUNCTION TRANS

      SUBROUTINE RESEED(NATOMS,P,RADIUS)
      ! {{{
      IMPLICIT NONE
      INTEGER J1, NATOMS
      DOUBLE PRECISION DPRAND, P(3*NATOMS), RADIUS, SR3, RANDOM
      
      SR3=DSQRT(3.0D0)
      DO J1=1,NATOMS
         RANDOM=(DPRAND()-0.5D0)*2.0D0
         P(3*(J1-1)+1)=RANDOM*DSQRT(RADIUS)/SR3
         RANDOM=(DPRAND()-0.5D0)*2.0D0
         P(3*(J1-1)+2)=RANDOM*DSQRT(RADIUS)/SR3
         RANDOM=(DPRAND()-0.5D0)*2.0D0
         P(3*(J1-1)+3)=RANDOM*DSQRT(RADIUS)/SR3
      ENDDO

      RETURN
      ! }}}
      ENDSUBROUTINE

!  The seed coordinates are at the end, not the beginning!!!!
      SUBROUTINE GSEED
      ! {{{
      USE COMMONS
      USE V

      IMPLICIT NONE
      
      DOUBLE PRECISION XMASS, YMASS, ZMASS, DIST, DUMMY, DMAX, &
     &                 DMIN, AX, BX
      INTEGER J1, I, J2


      OPEN(UNIT=10,FILE=SEED_FILE,STATUS='OLD')
      NSEED=0
      DO J1=1,NATOMS
         READ(10,*,END=10) COORDS(3*(NATOMS-J1)+1,1),COORDS(3*(NATOMS-J1)+2,1),COORDS(3*(NATOMS-J1)+3,1)
         NSEED=NSEED+1
      ENDDO
10    CLOSE(10)
      WRITE(LFH,'(A,I6,A)') 'Read core from file seed containing ',NSEED,' atoms'
      IF (FREEZECORE) THEN
         WRITE(LFH,'(A,I8,A)') 'Core will be fixed during the first ',NSSTOP,' quenches'
      ELSE
         WRITE(LFH,'(A,I8,A)') 'Core will be relaxed but reset for the first ',NSSTOP,' quenches'
      ENDIF
!
!  Centre the seed.
!
      IF (CENT) THEN
         XMASS=0.0D0
         YMASS=0.0D0
         ZMASS=0.0D0
         DO I=NATOMS,NATOMS-NSEED+1,-1
            XMASS=XMASS+COORDS(3*(I-1)+1,1)
            YMASS=YMASS+COORDS(3*(I-1)+2,1)
            ZMASS=ZMASS+COORDS(3*(I-1)+3,1)
         ENDDO
         XMASS=XMASS/NSEED
         YMASS=YMASS/NSEED
         ZMASS=ZMASS/NSEED
         DO I=NATOMS,NATOMS-NSEED+1,-1
            COORDS(3*(I-1)+1,1)=COORDS(3*(I-1)+1,1)-XMASS
            COORDS(3*(I-1)+2,1)=COORDS(3*(I-1)+2,1)-YMASS
            COORDS(3*(I-1)+3,1)=COORDS(3*(I-1)+3,1)-ZMASS
         ENDDO

!
!  Centre the other atoms.
!
         IF (NATOMS-NSEED.GT.1) THEN
            XMASS=0.0D0
            YMASS=0.0D0
            ZMASS=0.0D0
            DO I=1,NATOMS-NSEED
               XMASS=XMASS+COORDS(3*(I-1)+1,1)
               YMASS=YMASS+COORDS(3*(I-1)+2,1)
               ZMASS=ZMASS+COORDS(3*(I-1)+3,1)
            ENDDO
!           PRINT*,'NATOMS-NSEED=',NATOMS-NSEED
            XMASS=XMASS/(NATOMS-NSEED)
            YMASS=YMASS/(NATOMS-NSEED)
            ZMASS=ZMASS/(NATOMS-NSEED)
            DO I=1,NATOMS-NSEED
               COORDS(3*(I-1)+1,1)=COORDS(3*(I-1)+1,1)-XMASS
               COORDS(3*(I-1)+2,1)=COORDS(3*(I-1)+2,1)-YMASS
               COORDS(3*(I-1)+3,1)=COORDS(3*(I-1)+3,1)-ZMASS
            ENDDO
         ENDIF
      ENDIF
!
!  Find the largest radius vector of the seed.
!
      DIST=0.0D0
      DO J1=NATOMS,NATOMS-NSEED+1,-1
         DUMMY=COORDS(3*J1-2,1)**2+COORDS(3*J1-1,1)**2+COORDS(3*J1,1)**2
         IF (DUMMY.GT.DIST) DIST=DUMMY
      ENDDO
      DIST=DSQRT(DIST)
!
!  Shift the coordinates of the non-core atoms outside the core.
!
      DMAX=0.0D0
      DMIN=1.0D20
      DO J1=1,NATOMS-NSEED
         DUMMY=DSQRT(COORDS(3*J1-2,1)**2+COORDS(3*J1-1,1)**2+COORDS(3*J1,1)**2)
         IF (DUMMY.GT.DMAX) DMAX=DUMMY
         IF (DUMMY.LT.DMIN) DMIN=DUMMY
      ENDDO
      IF (DABS(DMAX-DMIN).GT.1.0D-5) THEN
         AX=1.2D0*DIST+DMIN*(1.2D0*DIST-DMIN-DMAX)/(DMAX-DMIN)
         BX=(DMIN+DMAX-1.2D0*DIST)/(DMAX-DMIN)
      ELSE
         AX=1.2D0*DIST
         BX=0.0D0
      ENDIF
      DO J1=1,NATOMS-NSEED
         DUMMY=DSQRT(COORDS(3*J1-2,1)**2+COORDS(3*J1-1,1)**2+COORDS(3*J1,1)**2)
         IF (DUMMY.LE.DIST) THEN
            COORDS(3*J1-2,1)=COORDS(3*J1-2,1)*(AX+BX*DUMMY)/DUMMY
            COORDS(3*J1-1,1)=COORDS(3*J1-1,1)*(AX+BX*DUMMY)/DUMMY
            COORDS(3*J1,1)  =COORDS(3*J1,1)  *(AX+BX*DUMMY)/DUMMY
         ENDIF
      ENDDO
      WRITE(LFH,75)
75    FORMAT('Coordinates:')
      WRITE(LFH,80) (COORDS(J1,1),J1=1,3*NATOMS)
80    FORMAT(3F15.5)
      IF (DUMPT) THEN
         WRITE(40,*) NATOMS
         WRITE(40,*) ' Initial coordinates'
         WRITE(40,45) (COORDS(J2,1),J2=1,3*(NATOMS-NSEED))
45       FORMAT('LA ',3F20.10)
         WRITE(40,46) (COORDS(J2,1),J2=3*(NATOMS-NSEED)+1,3*NATOMS)
46       FORMAT('LB',3F20.10)
      ENDIF

      DO J1=1,3*NATOMS
         COORDSO(J1,1)=COORDS(J1,1)
      ENDDO
      
      NS=NSEED

      RETURN
      ! }}}
      ENDSUBROUTINE
      
		FUNCTION PAIRDISTANCE(ATOM1,ATOM2)
		IMPLICIT NONE
		! PAIRDISTANCE is defined as it is the function name and hence the returned
		! value
		DOUBLE PRECISION :: PAIRDISTANCE 
		DOUBLE PRECISION, INTENT(IN) :: ATOM1(3),ATOM2(3)
		   
		PAIRDISTANCE=DSQRT((ATOM1(1)-ATOM2(1))**2+(ATOM1(2)-ATOM2(2))**2+(ATOM1(3)-ATOM2(3))**2)
		RETURN 
		END FUNCTION PAIRDISTANCE
        ! }}}
! GETRND SDPRND DPRAND {{{

! doxygen - GETRND {{{
!> @name         GETRND
! 
!> @brief        Get an array of random numbers inside the interval [XMIN,XMAX]
!
!> @param[in]    N              - dimension of the array RND
!> @param[out]   RND            - the generated array of random numbers 
!> @param[in]    XMIN,XMAX      
!
! }}}
!SUBROUTINE GETSEED(seed)
!integer seed


!ENDSUBROUTINE GETSEED

SUBROUTINE GETRND(RND,N,XMIN,XMAX)
! {{{
IMPLICIT NONE

! random number vector
DOUBLE PRECISION,DIMENSION(:),INTENT(OUT) :: RND
! range
DOUBLE PRECISION,INTENT(IN) :: XMIN,XMAX
DOUBLE PRECISION :: DX
! dimension of RND(:)
INTEGER,INTENT(IN) :: N
! seed
INTEGER :: ISEED
! loc
INTEGER I

DX=XMAX-XMIN

DO I=1,N 
        RND(I)=XMAX-DX*DPRAND()
ENDDO

RETURN
! }}}
END SUBROUTINE 

SUBROUTINE SDPRND (ISEED)
        ! declarations {{{
        DOUBLE PRECISION XMOD, YMOD, POLY(101), OTHER, OFFSET, X
        PARAMETER (XMOD = 1000009711.0D0, YMOD = 33554432.0D0)
        INTEGER ISEED, INDEX, IX, IY, IZ, I
        LOGICAL INITAL
        SAVE INITAL
        COMMON /RANDDP/ POLY, OTHER, OFFSET, INDEX
        DATA INITAL/.TRUE./
        ! }}}
        ! subroutine body {{{
!
!   ISEED should be set to an integer between 0 and 9999 inclusive;
!   a value of 0 will initialise the generator only if it has not
!   already been done.
!
        IF (INITAL .OR. ISEED .NE. 0) THEN
            INITAL = .FALSE.
        ELSE
            RETURN
        END IF
!
!   INDEX must be initialised to an integer between 1 and 101
!   inclusive, POLY(1...N) to integers between 0 and 1000009710
!   inclusive (not all 0), and OTHER to a non-negative proper fraction
!   with denominator 33554432.  It uses the Wichmann-Hill generator to
!   do this.
!
        IX = MOD(ABS(ISEED),10000)+1
        IY = 2*IX+1
        IZ = 3*IX+1
        DO 10 I = -10,101
            IF (I .GE. 1) POLY(I) = AINT(XMOD*X)
            IX = MOD(171*IX,30269)
            IY = MOD(172*IY,30307)
            IZ = MOD(170*IZ,30323)
            X = MOD(DBLE(IX)/30269.0D0+DBLE(IY)/30307.0D0+DBLE(IZ)/30323.0D0,1.0D0)
  10    CONTINUE
        OTHER = AINT(YMOD*X)/YMOD
        OFFSET = 1.0D0/YMOD
        INDEX = 1
        ! }}}
END SUBROUTINE SDPRND

FUNCTION DPRAND()
        ! DECLARATIONS {{{
        DOUBLE PRECISION XMOD, YMOD, XMOD2, XMOD4, TINY, POLY(101), DPRAND, &
         OTHER, OFFSET, X, Y
        PARAMETER (XMOD = 1000009711.0D0, YMOD = 33554432.0D0, &
        XMOD2 = 2000019422.0D0, XMOD4 = 4000038844.0D0,&
        TINY = 1.0D-17)
        INTEGER INDEX, N
        LOGICAL INITAL
        SAVE INITAL
        COMMON /RANDDP/ POLY, OTHER, OFFSET, INDEX
        DATA INITAL/.TRUE./
        ! }}}
        ! SUBROUTINE BODY {{{
!
!   THIS RETURNS A UNIFORM (0,1) RANDOM NUMBER, WITH EXTREMELY GOOD
!   UNIFORMITY PROPERTIES.  IT ASSUMES THAT DOUBLE PRECISION PROVIDES
!   AT LEAST 33 BITS OF ACCURACY, AND USES A POWER OF TWO BASE.
!
        IF (INITAL) THEN
            CALL SDPRND (0)
            INITAL = .FALSE.
        END IF
!
!   SEE [KNUTH] FOR WHY THIS IMPLEMENTS THE ALGORITHM DESCRIBED IN
!   THE PAPER.  NOTE THAT THIS CODE IS TUNED FOR MACHINES WITH FAST
!   DOUBLE PRECISION, BUT SLOW MULTIPLY AND DIVIDE; MANY, MANY OTHER
!   OPTIONS ARE POSSIBLE.
!
        N = INDEX-64
        IF (N .LE. 0) N = N+101
        X = POLY(INDEX)+POLY(INDEX)
        X = XMOD4-POLY(N)-POLY(N)-X-X-POLY(INDEX)
        IF (X .LT. 0.0D0) THEN
            IF (X .LT. -XMOD) X = X+XMOD2
            IF (X .LT. 0.0D0) X = X+XMOD
        ELSE
            IF (X .GE. XMOD2) THEN
                X = X-XMOD2
                IF (X .GE. XMOD) X = X-XMOD
            END IF
            IF (X .GE. XMOD) X = X-XMOD
        END IF
        POLY(INDEX) = X
        INDEX = INDEX+1
        IF (INDEX .GT. 101) INDEX = INDEX-101
!
!   ADD IN THE SECOND GENERATOR MODULO 1, AND FORCE TO BE NON-ZERO.
!   THE RESTRICTED RANGES LARGELY CANCEL THEMSELVES OUT.
!
   10   Y = 37.0D0*OTHER+OFFSET
        OTHER = Y-AINT(Y)
        IF (OTHER .EQ. 0.0D0) GO TO 10
        X = X/XMOD+OTHER
        IF (X .GE. 1.0D0) X = X-1.0D0
        DPRAND = X+TINY
        ! }}}
END FUNCTION DPRAND
! }}}
! INQF OPENF {{{
SUBROUTINE INQF(FILENAME,YESNO)
! {{{
LOGICAL,INTENT(OUT) :: YESNO
CHARACTER(LEN=*),INTENT(IN) :: FILENAME
CHARACTER(LEN=100) :: FLN

FLN=TRIM(ADJUSTL(FILENAME))
INQUIRE(FILE=FLN,EXIST=YESNO)
! }}}
END SUBROUTINE INQF

!> @name OPENF
!! @brief open files 
SUBROUTINE OPENF(FILEHANDLE,MODE,FILENAME)
! {{{

INTEGER, INTENT(IN) :: FILEHANDLE
CHARACTER (LEN=*), INTENT(IN) :: FILENAME
CHARACTER (LEN=*), INTENT(IN) :: MODE
CHARACTER (LEN=100) :: FLN

FLN=TRIM(ADJUSTL(FILENAME))

SELECTCASE(MODE)
        CASE(">")
                OPEN(UNIT=FILEHANDLE,FILE=FILENAME,STATUS="UNKNOWN",FORM="FORMATTED")
        CASE("O")
                OPEN(UNIT=FILEHANDLE,FILE=FILENAME,STATUS="OLD")
        CASE("<")
                OPEN(UNIT=FILEHANDLE,FILE=FILENAME,STATUS="OLD",ACTION="READ")
        CASE(">>")
                OPEN(FILEHANDLE,FILE=FILENAME,STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
        CASE("RW>>")
                OPEN(FILEHANDLE,FILE=FILENAME,STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND',ACTION="READWRITE")
        case("DA")
                OPEN(FILEHANDLE,FILE=FILENAME,ACCESS="DIRECT",STATUS='UNKNOWN',FORM='UNFORMATTED',RECL=8*3*NATOMS)
                
ENDSELECT
! }}}
ENDSUBROUTINE OPENF
! }}}
        ! gsort2 sort3 centre2 centrecom {{{
                ! gsort2 {{{ 
!     This subprogram performs a sort on the input data and
!     arranges it from smallest to biggest. The exchange-sort
!     algorithm is used.
!
      SUBROUTINE GSORT2(N,NATOMS)
      USE COMMONS, ONLY : QMINAV, QMINPCSMAV, CSMT
      USE V
      IMPLICIT NONE
      INTEGER NATOMS
      INTEGER J1, L, N, J3, J2, NTEMP
      DOUBLE PRECISION TEMP, C
!
      DO 20 J1=1,N-1
         L=J1
         DO 10 J2=J1+1,N
            IF (QMIN(L).GT.QMIN(J2)) L=J2
10       CONTINUE
         TEMP=QMIN(L)
         QMIN(L)=QMIN(J1)
         QMIN(J1)=TEMP
         IF (CSMT) THEN
            TEMP=QMINAV(L)
            QMINAV(L)=QMINAV(J1)
            QMINAV(J1)=TEMP
         ENDIF
         NTEMP=FF(L)
         FF(L)=FF(J1)
         FF(J1)=NTEMP
         DO J2=1,3*NATOMS
            C=QMINP(L,J2)
            QMINP(L,J2)=QMINP(J1,J2)
            QMINP(J1,J2)=C
         ENDDO
         IF (CSMT) THEN
            DO J2=1,3*NATOMS
               C=QMINPCSMAV(L,J2)
               QMINPCSMAV(L,J2)=QMINPCSMAV(J1,J2)
               QMINPCSMAV(J1,J2)=C
            ENDDO
         ENDIF
20    CONTINUE
      RETURN
      END SUBROUTINE GSORT2
      ! }}}
      ! sort3 {{{
!     This subprogram performs a sort on the input data and
!     arranges it from smallest to biggest. The exchange-sort
!     algorithm is used.
!
      SUBROUTINE SORT3(N,J3,A,B)
      IMPLICIT NONE
      INTEGER J1, L, N, J3, J2
      DOUBLE PRECISION TEMP, A(J3), B(3*J3)
!
      DO 20 J1=1,N-1
         L=J1
         DO 10 J2=J1+1,N
            IF (A(L).LT.A(J2)) L=J2
10       CONTINUE
         TEMP=A(L)
         A(L)=A(J1)
         A(J1)=TEMP
         DO J2=0,2
            TEMP=B(3*L-J2)
            B(3*L-J2)=B(3*J1-J2)
            B(3*J1-J2)=TEMP
         ENDDO
20    CONTINUE
      RETURN
      END SUBROUTINE
      ! }}}
      ! centre2 {{{
      SUBROUTINE CENTRE2(X)

      USE COMMONS
      USE V

      IMPLICIT NONE

      INTEGER I
      DOUBLE PRECISION,INTENT(INOUT),DIMENSION(NR) :: X
      DOUBLE PRECISION ::   RMASS(3)
      DOUBLE PRECISION ::   R(NATOMS,3)

      R=RESHAPE(X,(/ NATOMS,3 /))
      RMASS=SUM(R,DIM=1)/NATOMS
      do i=1,3
        R(1:NATOMS,i)=R(1:NATOMS,i)-RMASS(i)
      enddo
      X=PACK(R,.true.)
      RETURN
      END SUBROUTINE 
      ! }}}
      ! centrecom {{{
      SUBROUTINE CENTRECOM(X)

      USE COMMONS

      IMPLICIT NONE

      INTEGER I
      DOUBLE PRECISION X(3*NATOMS), TOTMASS,rmass(3)
      DOUBLE PRECISION ::   r(natoms,3)

      R=RESHAPE(X,(/ NATOMS,3 /))
      TOTMASS=SUM(MASSES)
      DO I=1,3
        RMASS(I)=SUM(R(1:NATOMS,I)*MASSES(1:NATOMS))/TOTMASS
      ENDDO
      do i=1,3
        R(1:NATOMS,i)=R(1:NATOMS,i)-RMASS(i)
      enddo
      X=PACK(R,.true.)
      
      IF (DEBUG) WRITE(LFH,'(A,3F15.10)') 'centre of mass reset to the origin from ',RMASS

      RETURN
      ENDSUBROUTINE CENTRECOM
      ! }}}
      ! }}}
      ! setvars {{{
SUBROUTINE SETVARS

! container radius
      IF (RADIUS.EQ.0.0D0) THEN
         RADIUS=2.0D0+(3.0D0*NATOMS/17.77153175D0)**(1.0D0/3.0D0)
         IF (P46) THEN
            RADIUS=RADIUS*3.0D0
         ELSE 
            RADIUS=RADIUS*2.0D0**(1.0D0/6.0D0)
         ENDIF
      ENDIF

      RADIUS=RADIUS**2

END SUBROUTINE SETVARS
! }}}
! initvars {{{
RECURSIVE SUBROUTINE INITVARS(S)

INTEGER NPCALL
COMMON /PCALL/ NPCALL
character(len=*) S
! subroutine body {{{

selectcase(S)
    case("LOG")
! logicals {{{
LBFGST=.TRUE.
USESUF=.FALSE.
PULLT=.TRUE.
P46=.FALSE.
G46=.TRUE.
BLNT=.FALSE.
MYBLNT=.FALSE.
TARGET=.FALSE.
TRACKDATAT=.FALSE.
DEBUG=.FALSE.
DUMPT=.FALSE.
NORESET=.FALSE.
RMST=.FALSE.
      
SAVEQ=.TRUE.
SORTT=.FALSE.

ALLOCATE(FIXSTEP(1),FIXTEMP(1),FIXBOTH(1))
FIXSTEP=.FALSE.
FIXTEMP=.FALSE.
FIXBOTH=.FALSE.
! whether we are doing a final quench
!FQFLAG=.FALSE.
! }}}
    case("FILES")
! files {{{
IFH=50
LFH=IFH
EA_FH=IFH+1
ENERGY_FH=IFH+2
MARKOV_FH=IFH+3
BEST_FH=IFH+4
PAIRDIST_FH=IFH+5
COORDS_FH=IFH+6
DATA_FH=IFH+7
RMSD_FH=IFH+8
LE_FH=IFH+9

C_FILE="coords"
D_FILE="data"
LE_FILE="le"
E_FILE="e.tex"
O_FILE="out"
SEED_FILE="seed"
EA_FILE="ea"

IF (USESUF) THEN
	LE_FILE=adjustr(SUF)//"le"
	E_FILE=SUF//"e.tex"
	EA_FILE=SUF//"ea"
	O_FILE=adjustr(SUF)//"out"
ENDIF

! }}}
    case("VARS")
! other {{{

      BLNTYPE="GO"

      IF (P46) THEN
        BLNTYPE="WT"
      ELSEIF(G46)THEN
        BLNTYPE="GO"
      ENDIF
      NRG=1
            
TFAC=1.0D0          ! temperature multiplier

NPCALL=0
NSEED=0
NS=0

NSAVE=10            ! number of saved lowest-energy geometries
NATOMS=46           ! number of particles
!EAMP=0.01D0         ! exponential multiplier


PATOM1=1            ! atom 1 for pulling 
PATOM2=NATOMS       ! atom 2 for pulling

MUPDATE=4           ! Number of LBFGS updates
MAXBFGS=0.4D0       ! Maximum BFGS step size
DGUESS=0.1D0        ! DGUESS: Guess for initial diagonal elements in LBFGS

FQMAX=1.0D-5        ! FQMAX: same meaning as for SQMAX, but for final quenches only.
SQMAX=1.0D-3        ! SQMAX: convergence criterion for the RMS force in the basin-hopping quenches.
                    ! note: used in QUENCH() 
ECONV=0.02D0
RADIUS=0.0D0
! maximum number of iterations allowed in conjugate gradient searches
MAXIT=500
! Maximum allowed energy rise during a minimisation
MAXERISE=1.0D-10
! Maximum allowed energy fall during a minimisation
MAXEFALL=-HUGE(ONE)
! Used in ACCREJ
!FAC0=1.05D0
! "Fixing" option (regarding STEP, TEMP and accept ratio for quenches) 
!FIXOPT='T'

ARMA=0.4D0
ARMB=0.4D0
EPSSPHERE=0.0D0
          
NACCEPT=50
NRELAX=0
! }}}
    case("ARR")
! {{{
    ALLOCATE(ACCRAT(1),STEP(1),OSTEP(1),ASTEP(1),TEMP(1))
    ALLOCATE(BEADLETTER(NATOMS)) 
    ACCRAT=0.5D0
    TEMP=0.035D0        ! Temperature
	STEP=0.3D0
	OSTEP=0.3D0
	ASTEP=0.3D0
	MCSTEPS=10000
      ! hydrophobic - B 1
      ! neutral N    3
      ! hydrophilic L 2
      BEADLETTER(1:9)="B"
      BEADLETTER(10:12)="N"
      BEADLETTER(13:19:2)="L"
      BEADLETTER(14:20:2)="B"
      BEADLETTER(21:23)="N"
      BEADLETTER(24:32)="B"
      BEADLETTER(33:35)="N"
      BEADLETTER(36:46:2)="L"
      BEADLETTER(37:45:2)="B"
! }}}
    case("ALL")
      CALL INITVARS("LOG")
      CALL INITVARS("VARS")
      CALL INITVARS("ARR")
      CALL SETVARS
endselect
! }}}
END SUBROUTINE INITVARS
! }}}
! r_coords {{{
! doxygen R_COORDS {{{
!> @name RCOORDS
!> @brief Generate a random set of coordinates, based on:
!> @param[in] NATOMS number of particles
!> @param[in] RADIUS container radius
!> @param[out] COORDS  randomly generated coordinates
! }}}
SUBROUTINE R_COORDS(NATOMS,RADIUS,COORDS,ISEED)
! declarations {{{
! subroutine parameters 
INTEGER, INTENT(IN) :: NATOMS
DOUBLE PRECISION,INTENT(IN) :: RADIUS 
DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: COORDS
DOUBLE PRECISION,DIMENSION(3*NATOMS) :: RND
INTEGER,INTENT(IN) :: ISEED

! local parameters 
INTEGER NR
DOUBLE PRECISION :: SR3
! }}}
! {{{
NR=3*NATOMS
SR3=DSQRT(3.0D0)
CALL GETRND(RND,NR,-1.0D0,1.0D0)
COORDS=RADIUS*RND/SR3
! }}}
END SUBROUTINE R_COORDS
! }}}

      SUBROUTINE COUNTATOMS
!op226> Declarations {{{ 
      USE COMMONS, ONLY : NATOMS

      IMPLICIT NONE

      INTEGER :: EOF
      LOGICAL :: YESNO
!op226>}}} 
      ! {{{

! commented  {{{
!
!  If the current working directory contains more than one of these files
!  then the precedence is coords, then input.crd, then coords.amber
!  OPTIM does this a bit better by calling getparams first to see if
!  we are actually doing AMBER or CHARMM. 
!

      !YESNOA=.FALSE.
      !YESNOAMH=.FALSE.
      !YESNOA9=.FALSE.
      !INQUIRE(FILE='pro.list',EXIST=YESNOAMH)
      !INQUIRE(FILE='coords.amber',EXIST=YESNOA)
      !INQUIRE(FILE='input.crd',EXIST=YESNOC)
      !INQUIRE(FILE='coords.inpcrd',EXIST=YESNOA9)
      ! }}}

      YESNO=.FALSE.
      INQUIRE(FILE=C_FILE,EXIST=YESNO)

      IF (YESNO) THEN
         OPEN(UNIT=COORDS_FH,FILE=C_FILE,STATUS='OLD')
         NATOMS=0
         REWIND(COORDS_FH)
         DO
            READ(COORDS_FH,*,IOSTAT=EOF)
            IF (EOF==0) THEN
               NATOMS = NATOMS + 1
            ELSE
               EXIT
            ENDIF
         ENDDO
         CLOSE(COORDS_FH)
      ELSE
         PRINT '(A)','ERROR - no coords, input.crd, coords.inpcrd or coords.amber file'
         STOP
      ENDIF

      NR=3*NATOMS

      CLOSE(7)
      ! }}}
      END SUBROUTINE COUNTATOMS

      ! getmodel printvars printhelp {{{
SUBROUTINE GETMODEL
! {{{

IF (P46) THEN
  MODEL="P46 THREE-COLOUR OFF-LATTICE PROTEIN MODEL, WILD-TYPE"
ELSEIF(G46) THEN
  MODEL="P46 THREE-COLOUR OFF-LATTICE PROTEIN MODEL, GO-LIKE"
ELSEIF(BLNT) THEN
  MODEL="GENERAL BLN MODEL"
ENDIF
! }}}
ENDSUBROUTINE GETMODEL

!! print variables 
SUBROUTINE PRINTVARS(FH)
! {{{
! vars                                                                       {{{
! sub
integer fh
! loc
character(100) fmt(10)
character(40) s
integer i
! }}}
include '../include/fmt.i.f90'

CALL ED(fh)
write(fh,10) "PARAMETER VALUES" !                                             {{{

CALL GETMODEL
WRITE(fh,*) "MODEL" 
write(fh,*) TRIM(MODEL)
CALL ED(fh)
write(fh,1)  "PARAMETER DESCRIPTION",         "NAME",                  "VALUE"
! pd:general                                                                 {{{
write(fh,11) s_stars
write(fh,10) "GENERAL" 
write(fh,11) s_stars
!s=adjustl(s)
write(fh,3)  "Number of particles",                          "NATOMS",     NATOMS
write(fh,2)  "Container radius",                             "RADIUS",     RADIUS
write(fh,3)  "Number of saved lowest energy geometries",     "NSAVE",      NSAVE
write(fh,3)  "Number of basin-hopping steps",                "MCSTEPS",    MCSTEPS
write(fh,2)  "Temperature",                                  "TEMP",       TEMP
write(fh,2)  "Acceptance ratio",                             "ACCRAT",     ACCRAT
write(fh,2)  "Energy difference criterion for minima",       "ECONV",      ECONV
write(fh,2)  "Final quench tolerance for RMS gradient ",     "FQMAX",      FQMAX
write(fh,2)  "Quench convergence criterion for RMS gradient ", "SQMAX",    SQMAX
write(fh,3)  "Maximum number of iterations",   "MAXIT", MAXIT
write(fh,11) "(sloppy and final quenches)"
write(fh,2)  "Temperature multiplier",                       "TFAC",       TFAC
write(fh,3)  "",                                             "NACCEPT",    NACCEPT
write(fh,3)  "",                                             "NRELAX",     NRELAX
write(fh,11) s_stars !                                                        }}}
! pd:lbfgs                                                                     {{{
write(fh,10) "LBFGS parameters"
write(fh,11) s_stars
write(fh,3) "Number of LBFGS updates",   "MUPDATE",  MUPDATE
write(fh,2) "Maximum BFGS step size",   "MAXBFGS",   MAXBFGS
write(fh,2) "Guess for initial diagonal elements in BFGS", "DGUESS", DGUESS
! }}} 
! }}}
CALL ED(fh)
! }}}
END SUBROUTINE PRINTVARS

SUBROUTINE PRINTHELP(FH)
INTEGER FH

CALL ED(FH)
write(fh,*) ''
write(fh,*) 'Command-line options:'
write(fh,*) ''
write(fh,*) '   -ca NACCEPT'
write(fh,*) '   -f  PFORCE'
write(fh,*) '   -g  turn on debug printing'
write(fh,*) ''
write(fh,*) '   -v  display version'
write(fh,*) '   -pv display default parameter values'
write(fh,*) ''
write(fh,*) '   -h  print this help message'
write(fh,*) '   -s SUF Output files with suffix SUF, e.g., out => SUF.out'
write(fh,*) '   -steps STEPS Number of Monte-Carlo steps'
write(fh,*) ''
CALL ED(FH)

ENDSUBROUTINE PRINTHELP 
! }}}
    ! am {{{
    SUBROUTINE AM(S)
        CHARACTER(LEN=*) S

        selectcase(S)
            CASE("MAIN")
               ALLOCATE(MSCREENC(3*NATOMS),VT(NATOMS))
               ALLOCATE(FF(NSAVE),QMIN(NSAVE))
               ALLOCATE(QMINP(NSAVE,3*NATOMS))
               ALLOCATE(COORDSO(3*NATOMS,1),VAT(NATOMS,1),VATO(NATOMS,1),COORDS(3*NATOMS,1))
            CASE("INIT")
      ALLOCATE(NQ(1),EPREV(1))
      ALLOCATE(NCORE(1))
        endselect
    ENDSUBROUTINE AM
    ! }}}


      END MODULE F
