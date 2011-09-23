      MODULE F

      USE V 

      IMPLICIT NONE 

      CONTAINS

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

      ! interfaces 
      INTERFACE
        SUBROUTINE GSAVEIT(EREAL,P,NP)
            INTEGER NP
            DOUBLE PRECISION :: EREAL
            DOUBLE PRECISION,DIMENSION(:) ::   P
        END SUBROUTINE GSAVEIT

        SUBROUTINE TRANSITION(ENEW,EOLD,ATEST,NP,RANDOM,MCTEMP)
            DOUBLE PRECISION ENEW, EOLD, RANDOM, MCTEMP
            LOGICAL ATEST
            INTEGER NP
        END SUBROUTINE TRANSITION

        SUBROUTINE ACCREJ(NSUCCESS,NFAIL,JP,NSUCCESST,NFAILT)
          INTEGER NSUCCESS(NPAR), NFAIL(NPAR), JP, NFAILT(NPAR), NSUCCESST(NPAR)
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

      END MODULE F
