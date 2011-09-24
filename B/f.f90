      MODULE F

      USE V 
      USE COMMONS

      IMPLICIT NONE 

      ! interfaces  {{{
      INTERFACE

        DOUBLE PRECISION FUNCTION DPRAND()
        END FUNCTION DPRAND

        SUBROUTINE MYRESET(JP,NATOMS,NPAR,NSEED)
            INTEGER JP,NATOMS,NPAR,NSEED
        END SUBROUTINE MYRESET

        SUBROUTINE GSAVEIT(EREAL,P,NP)
            INTEGER NP
            DOUBLE PRECISION :: EREAL
            DOUBLE PRECISION,DIMENSION(:) ::   P
        END SUBROUTINE GSAVEIT

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

        SUBROUTINE MYLBFGS(N,M,XCOORDS,DIAGCO,EPS,MFLAG,ENERGY,ITMAX,ITDONE,RESET,NP)

	      INTEGER :: N,M,ITMAX,ITDONE,NP
	      DOUBLE PRECISION,DIMENSION(:) ::   XCOORDS
	      LOGICAL DIAGCO,MFLAG,RESET
	      DOUBLE PRECISION ::   EPS, ENERGY

        END SUBROUTINE MYLBFGS

        SUBROUTINE QUENCH(QTEST,NP,ITER,TIME,BRUN,QDONE,P)

	      DOUBLE PRECISION, DIMENSION(:) :: P   
	      INTEGER NP,ITER,BRUN,QDONE
	      LOGICAL QTEST
	      DOUBLE PRECISION ::   TIME

        END SUBROUTINE QUENCH

      SUBROUTINE POTENTIAL(X,GRAD,EREAL,GRADT,SECT)
	      DOUBLE PRECISION, DIMENSION(:),INTENT(IN) :: X
	      DOUBLE PRECISION, DIMENSION(:),INTENT(OUT) :: GRAD
	      DOUBLE PRECISION, INTENT(OUT) :: EREAL
	      LOGICAL, INTENT(IN) :: GRADT, SECT       
      END SUBROUTINE POTENTIAL

      END INTERFACE
      ! }}}

      CONTAINS

      ! trans gseed reseed pairdistance  {{{

      DOUBLE PRECISION FUNCTION TRANS(X,XMIN,GAMMA)
      IMPLICIT NONE
      DOUBLE PRECISION X, XMIN, GAMMA

      TRANS=1.0D0 - EXP(-GAMMA*(X-XMIN))

      RETURN
      END 

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
      END

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
      END
      
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

      SUBROUTINE SETVARS
!{{{
      USE COMMONS, ONLY : P46,G46

      BLNTYPE="GO"

      IF (P46) THEN
        BLNTYPE="WT"
      ELSEIF(G46)THEN
        BLNTYPE="GO"
      ENDIF
!}}}
      END SUBROUTINE SETVARS

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

         OPEN(UNIT=7,FILE=C_FILE,STATUS='OLD')
         DO
            READ(7,*,IOSTAT=EOF)
            IF (EOF==0) THEN
               NATOMS = NATOMS + 1
            ELSE
               EXIT
            ENDIF
         ENDDO
        !ELSEIF (YESNOAMH) THEN 
!        ! {{{
         !open(unit=30,file='pro.list',status='old',form='formatted')
         !read (30,1000)tarfl
!1000     format(a5)
         !close(30)

          !open(30,file='proteins/'//tarfl,status='old')
            !read(30,*)
            !read(30,*)nres
            !if (nres.gt.500) then
                !write(6,*) 'failure nres gr than 500 countatoms'
                !stop
            !endif
            !read (30,25)(seq(i_res),i_res=1,nres)
!!            write(6,25)(seq(i_res),i_res=1,nres)
!25         format(25(i2,1x))
          !close(30)

          !NOGLY = 0
          !GLY = 0

           !do i_res=1,nres
             !if (seq(i_res).ne.8) NOGLY = NOGLY +1
             !if (seq(i_res).eq.8) GLY = GLY +1
           !enddo

            !Number_of_Atoms = NOGLY*3 + GLY*2
      !ELSE IF (YESNOA9) THEN
!!         OPEN(UNIT=7,FILE='coords.gayberne',STATUS='OLD')
!!         PRINT '(A)','reading coordinates from file coords.gayberne'

         !inpcrd1='coords.inpcrd'
!!         inpcrd1=trim(adjustl(inpcrd1))
         !call amberinterface(Number_of_Atoms,1,inpcrd1,LFH)

      !ELSEIF (YESNOC) THEN
         !OPEN(UNIT=7,FILE='input.crd',STATUS='OLD')
         !do
           !read(7,*) myline
           !if (myline(1:1)=='*') then ! SAT This is the goddamn CHARMM comment line
              !cycle
           !else
              !read(myline,*) Number_of_Atoms
              !exit
           !endif
         !enddo

!! DAE We also need to find out what MAXAIM is in CHARMM, and set MXATMS in OPTIM to be the same, so that those arrays which
!! are passed between the two can be declared correctly. MXATMS is now stored in modmxatms.

         !CALL GETMAXAIM
         !WRITE(LFH,'(A,I8)') 'countatoms> Number_of_Atoms=',Number_of_Atoms
      !ELSEIF (YESNOA) THEN
         !OPEN(UNIT=7,FILE='coords.amber',STATUS='OLD')
         !do
            !read(7,'(A3)',iostat=eof) check
            !if (eof.LT.0) then
               !PRINT *,'End of file before all information specified'
               !STOP
            !ENDIF
            !IF (check.EQ.'end' .OR. check.EQ.'END' .OR. check.EQ.'End') THEN
               !CLOSE(7)
               !EXIT
            !ENDIF
            !Number_of_Atoms = Number_of_Atoms + 1
         !enddo
         !! }}}
      ELSE
         PRINT '(A)','ERROR - no coords, input.crd, coords.inpcrd or coords.amber file'
         STOP
      ENDIF

      CLOSE(7)
      ! }}}
      END SUBROUTINE COUNTATOMS

      END MODULE F
