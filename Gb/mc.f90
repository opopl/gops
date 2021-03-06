      MODULE MCFUNC

      USE COMMONS 
      USE QMODULE
      USE FUNC

      IMPLICIT NONE

      CONTAINS

     ! REST {{{
      SUBROUTINE REST(ITERATIONS,TIME,J1,RCOORDS,RMIN,RVAT,JACCPREV)
        ! dec {{{
      USE COMMONS

      IMPLICIT NONE

      INTEGER ITERATIONS, J2, JACCPREV, J1, NQTOT, BRUN, QDONE
      DOUBLE PRECISION TIME, POTEL, RCOORDS(3*NATOMS), RMIN, RVAT(NATOMS), SCREENC(3*NATOMS)
      COMMON /MYPOT/ POTEL
      COMMON /TOT/ NQTOT
      ! }}}
      ! body {{{

!10    CALL HSMOVE(COORDS(1:3*NATOMS,1:NPAR),1,NHSRESTART)
!  next line should be uncommented if routine is made availabe to use with CHARMM
      CALL QUENCH(.FALSE.,1,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
!
!  Bad idea to accept this quench configuration unconditionally - it could be unconvergeable.
!
      IF (POTEL-EPREV.GT.10.0D0*ABS(EPREV)) THEN
         DO J2=1,3*NATOMS
            COORDS(J2)=COORDSO(J2)
         ENDDO
         !GOTO 10
      ENDIF
      JACCPREV=J1
      NQTOT=NQTOT+1
      WRITE(LFH,'(A,I6,A)') ' Restarting using ',NHSRESTART,' hard sphere moves'
      WRITE(LFH,'(A,I7,A,F20.10,A,I5,A,G12.5,A,F20.10,A,F11.1)') 'Restart Qu ',NQ,' E=',&
     &              POTEL,' steps=',ITERATIONS,' RMS=',RMS,' t=',TIME-TSTART
      DO J2=1,3*NATOMS
         COORDSO(J2)=COORDS(J2)
         RCOORDS(J2)=COORDS(J2)
      ENDDO
      DO J2=1,NATOMS
         VATO(J2)=VAT(J2)
         RVAT(J2)=VAT(J2)
      ENDDO
      EPREV=POTEL
      RMIN=POTEL

      RETURN
      ! }}}
      END SUBROUTINE
     !     }}}

      SUBROUTINE MC(NSTEPS,SCALEFAC,SCREENC)
!OP226> DECLARATIONS {{{ 

      USE COMMONS
      USE QMODULE , ONLY : QMIN, QMINP, INTEQMIN
      USE PORFUNCS

      IMPLICIT NONE

      ! SUBROUTINE 
      INTEGER ::    NSTEPS
      DOUBLE PRECISION ::   SCALEFAC,SCREENC(3*NATOMS)

      ! LOCAL {{{
      INTEGER ::    ISTEP
      INTEGER :: NSUCCESS, NFAIL, NFAILT, NSUCCESST
      INTEGER :: ITERATIONS
      INTEGER :: NQTOT,JACCPREV

      ! USED IN QUENCH, NEWRES
      INTEGER BRUN 
      INTEGER QDONE,NDONE
      INTEGER J1, J2, JP, J5, I, J, J3, J4
      INTEGER :: JBEST
      DOUBLE PRECISION :: EBEST,EPPREV
      DOUBLE PRECISION,DIMENSION(3*NATOMS) :: GRAD,SAVECOORDS,TEMPCOORDS,RCOORDS,RCOORDSO
      DOUBLE PRECISION,DIMENSION(3*NATOMS) :: BESTCOORDS
      DOUBLE PRECISION,DIMENSION(NATOMS) :: RVAT,RVATO
      DOUBLE PRECISION ::   EPSSAVE
      DOUBLE PRECISION :: POTEL, RANDOM, DPRAND
      DOUBLE PRECISION ::  TIME   
      DOUBLE PRECISION ::  RMIN, RMINO
      DOUBLE PRECISION ::  T,  E, ER, W 
      DOUBLE PRECISION :: OPOTEL  
!  CSW34> PAIRDIST VARIABLES
      INTEGER :: PAIRCOUNTER
      DOUBLE PRECISION :: ATOM1(3),ATOM2(3)

      LOGICAL EVAP, ATEST, STAY, EVAPREJECT, LOPEN
      COMMON /EV/ EVAP, EVAPREJECT
      COMMON /MYPOT/ POTEL
      COMMON /TOT/ NQTOT
!OP226>}}} 
      ! }}}
!OP226> SUBROUTINE BODY {{{ 

      ! intro{{{
      INQUIRE(UNIT=1,OPENED=LOPEN)

      IF (LOPEN) THEN
         WRITE(*,'(A,I2,A)') 'MC> A ERROR *** UNIT ', 1, ' IS NOT FREE '
         STOP
      ENDIF

      EVAPREJECT=.FALSE.

      INQUIRE(UNIT=1,OPENED=LOPEN)

      IF (LOPEN) THEN
         WRITE(*,'(A,I2,A)') 'MC> B ERROR *** UNIT ', 1, ' IS NOT FREE '
         STOP
      ENDIF

      NDONE=0
      NQ=NDONE

      IF (NACCEPT.EQ.0) NACCEPT=NSTEPS+1

      JACCPREV=0
      NQTOT=0
      RMINO=1.0D100
      RMIN=1.0D100

         TMOVE=.TRUE.
         OMOVE=.TRUE.
         NSUCCESS=0
         NFAIL=0
         NSUCCESST=0
         NFAILT=0
      ! }}}

!  CALCULATE THE INITIAL ENERGY AND SAVE IN EPREV
!OP226>{{{ 
      WRITE(LFH,'(A)') 'CALCULATING INITIAL ENERGY'
      EPSSAVE=EPSSPHERE
      EPSSPHERE=0.0D0
      DO JP=1,NPAR
         CALL QUENCH(.FALSE.,JP,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
         NQTOT=NQTOT+1
         WRITE(LFH,'(A,I10,A,F20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)') &
            & 'QU ',NQ,&
            & ' E=',POTEL,&
            & ' STEPS=',ITERATIONS,&
            & ' RMS=',RMS,&
            & ' MARKOV E=',POTEL,&
            & ' T=',TIME-TSTART

!  EPREV SAVES THE PREVIOUS ENERGY IN THE MARKOV CHAIN.
!  EBEST AND JBEST RECORD THE LOWEST ENERGY SINCE THE LAST RESEEDING AND THE
!  STEP IT WAS ATTAINED AT. BESTCOORDS CONTAINS THE CORRESPONDING COORDINATES.
 
         EPREV=POTEL
         EPPREV=0.0D0
         EBEST=POTEL
         BESTCOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS)
         JBEST=0
         RMIN=POTEL
         RCOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS)
         COORDSO(1:3*NATOMS)=COORDS(1:3*NATOMS)
         VATO(1:NATOMS)=VAT(1:NATOMS)
      ENDDO

      EPSSPHERE=EPSSAVE
!OP226>}}} 

!OP226> GMIN_OUT: STARTING MC RUN ...; TEMPERATURE WILL ... {{{ 
      WRITE(LFH,'(A,I10,A)') 'STARTING MC RUN OF ',NSTEPS,' STEPS'
      WRITE(LFH,'(A,F15.8,A)') 'TEMPERATURE WILL BE MULTIPLIED BY ',SCALEFAC,' AT EVERY STEP'
!OP226>}}} 


!  MAIN BASIN-HOPPING LOOP 
! {{{
      DO J1=NDONE+1,NSTEPS 
         ISTEP = J1

         CALL FLUSH(LFH)
!
!  ********************************* LOOP OVER NPAR PARALLEL RUNS ******************************
!
         DO JP=1,NPAR 
!  MAM1000> THE DEFAULT TEMPERATURE USED FOR THE MC ACCEPTANCE CRITERION IS THE ONE DERIVED IN
!           THE INITIALISATION SECTION ABOVE.  MCTEMP IS PASSED TO THE SUBROUTINE TRANSITION, WHERE
!           THE ACCEPTANCE/REJECTION DECISION IS MADE.  HOWEVER, INDIVIDUAL MC MOVES CAN OVERRIDE
!           THIS TEMPERATURE BY SETTING MCTEMP TO A DIFFERENT VALUE FOR THE CURRENT STEP.
            MCTEMP = TEMP
!  ORDINARY STEPS.
 
23          CONTINUE
!
! DON'T CALL SYMMETRY UNLESS THE MINIMUM IN THE MARKOV CHAIN HAS CHANGED.
! WE SHOULD REALLY CHECK IF THE MINIMUM HAS CHANGED SINCE THE LAST CALL TO SYMMETRY,
! WHICH CAN BE DONE WITH ABS(ELASTSYM(JP)-EPREV) IF NSYMINTERVAL=1.
!
! CSW34> COORDINATES ARE SAVED SO THAT MOVES CAN BE UNDONE
           SAVECOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS)
           CALL TAKESTEP
           NQ=NQ+1
           CALL QUENCH(.FALSE.,JP,ITERATIONS,TIME,BRUN,QDONE,SCREENC)  
           NQTOT=NQTOT+1
!  OUTPUT
           WRITE(LFH,'(A,I10,A,F20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)') 'QU ',NQ,' E=',&
     &                 POTEL,' STEPS=',ITERATIONS,' RMS=',RMS,' MARKOV E=',EPREV,' T=',TIME-TSTART
          CALL FLUSH(LFH)

          ! PAIRDIST TRACKDATA {{{
          IF (PAIRDISTT) THEN
            ! {{{
!     WRITE END OF PREVIOUS LINE AS USING ADVANCE="NO"
                WRITE(MYPUNIT,*) " "
!     WRITE CURRENT QUENCH NUMBER
                WRITE(MYPUNIT,'(I10)',ADVANCE="NO") NQ
!     FOR EACH PAIR, ASSIGN ATOM1 AND ATOM2 ARRAYS CONTAINING COORDINATES
             DO PAIRCOUNTER=1,NPAIRS
                ATOM1(:)=COORDS(3*PAIRDIST(PAIRCOUNTER,1)-2:3*PAIRDIST(PAIRCOUNTER,1))
                ATOM2(:)=COORDS(3*PAIRDIST(PAIRCOUNTER,2)-2:3*PAIRDIST(PAIRCOUNTER,2))
!     CALL PAIRDISTANCE WITH (X,Y,Z) FOR EACH ATOM
                WRITE(MYPUNIT,'(F10.4)',ADVANCE="NO") PAIRDISTANCE(ATOM1,ATOM2) 
             ENDDO
             FLUSH(MYPUNIT)
             ! }}}
          ENDIF
!
          IF (TRACKDATAT) THEN
            ! {{{
             WRITE(MYEUNIT,'(I10,F20.10)') J1,POTEL
             WRITE(MYMUNIT,'(I10,G20.10)') J1,EPREV
             WRITE(MYBUNIT,'(I10,G20.10)') J1,QMIN(1)
             CALL FLUSH(MYEUNIT)
             CALL FLUSH(MYMUNIT)
             CALL FLUSH(MYBUNIT)
             ! }}}
          ENDIF
          ! }}}

!  CHECK FOR RESEEDING.

            IF (EVAP .AND. .NOT.EVAPREJECT) THEN
               NFAIL=NFAIL+1
               CALL MYRESET(NATOMS,NSEED)
               IF (DEBUG) THEN
                  WRITE(LFH,33) JP,J1,POTEL,EPREV,NSUCCESS,NFAIL
33                FORMAT('JP,J1,POTEL,EPREV,NSUC,NFAIL=',I2,I6,2F15.7,2I6,' EVAP,REJ')
               ENDIF
            ELSE
!     CSW34> A SERIES OF TESTS START HERE TO CHECK IF A STRUCTURE SHOULD
!     BE ALLOWED INTO THE MARKOV CHAIN. IF IT FAILS A TEST, THE ATEST
!     VARIABLE WILL BE SET TO .FALSE. 
               ATEST=.TRUE. 

               IF (ATEST) THEN
                  CALL TRANSITION(POTEL,EPREV,ATEST,JP,RANDOM,MCTEMP)
               ENDIF

!  SANITY CHECK TO MAKE SURE THE MARKOV ENERGY AGREES WITH COORDSO. 
!  STOP IF NOT TRUE.

               IF (DEBUG.OR.CHECKMARKOVT) THEN 
                  ! {{{
                  CALL POTENTIAL(COORDSO(:),GRAD,OPOTEL,.FALSE.,.FALSE.)
                  IF (ABS(OPOTEL-EPREV).GT.EDIFF) THEN
                     IF (EVAP) THEN
                        WRITE(LFH,'(3(A,G20.10))') &
                        & 'MC> WARNING - ENERGY FOR SAVED COORDINATES ',OPOTEL,&
                        & ' DIFFERS FROM MARKOV ENERGY ',EPREV,&
                        & ' BECAUSE AN ATOM MOVED OUTSIDE THE CONTAINER'
                     ELSE
                        WRITE(LFH,'(2(A,G20.10))') &
                        'MC> ERROR - ENERGY FOR COORDINATES IN COORDSO=',OPOTEL,&
                        & ' BUT MARKOV ENERGY=',EPREV
                        STOP
                     ENDIF
                  ENDIF
                     ! }}}
               ENDIF

! ACCEPT OR REJECT STEP. IF THE QUENCH DID NOT CONVERGE THEN ALLOW A
! POTENIAL MOVE, BUT COUNT IT AS A REJECTION IN TERMS OF NSUCCESS AND
! NFAIL. THIS WAY WE WILL ACCEPT A LOWER MINIMUM IF FOUND, BUT THE STEPS WON;T BECOME SO BIG.
! HOWEVER, FOR TIP5P SOME COLD FUSION EVENTS THAT HAD NOT ACTUALLY REACHED THE THRESHOLD FOR
! REJECTION WERE ACTUALLY ACCEPTED. MUST PREVENT THIS!
               IF (ATEST) THEN
                 ! {{{

                  IF (DEBUG) THEN
                     WRITE(LFH,34) JP,RANDOM,POTEL,EPREV,NSUCCESS,NFAIL
34                   FORMAT('JP,RAN,POTEL,EPREV,NSUC,NFAIL=',I2,3F15.7,2I6,' ACC')
                  ENDIF

                  IF ((J1-JACCPREV.GT.NRELAX).AND.ABS(POTEL-EPREV).GT.EDIFF) THEN
!                    NRELAX=J1-JACCPREV
!                    IF (RESTART) WRITE(LFH,'(A,I6,A)') ' RELAXATION TIME SET TO ',NRELAX,' STEPS'
                     JACCPREV=J1
                  ENDIF

                  IF (QDONE.EQ.1) THEN
                     NSUCCESS=NSUCCESS+1
                  ELSE
                     NFAIL=NFAIL+1
                  ENDIF

                  EPPREV=EPREV
                  EPREV=POTEL
                  COORDSO(1:3*NATOMS)=COORDS(1:3*NATOMS)
                  VATO(1:NATOMS)=VAT(1:NATOMS)
                  ! }}}
               ELSE
                  ! {{{
                  NFAIL=NFAIL+1
                  CALL MYRESET(NATOMS,NSEED)
                  IF (DEBUG) THEN
                     WRITE(LFH,36) RANDOM,POTEL,EPREV,NSUCCESS,NFAIL
36                   FORMAT('RAN,POTEL,EPREV,NSUC,NFAIL=',I2,3F15.7,2I6,' REJ')
                  ENDIF
                  ! }}}
               ENDIF
            ENDIF
!
!  IF RESTART THEN RESEED IF WE HAVEN T ACCEPTED A STEP IN TWICE THE RELAXATION TIME.
!
            IF (RESTART.AND.(J2-JACCPREV.GT.1.1D0*NRELAX)) &
                & CALL REST(ITERATIONS,TIME,J1,RCOORDS,RMIN,RVAT,JACCPREV)
!
!  CHECK THE ACCEPTANCE RATIO.
! 
            IF ((MOD(J1,NACCEPT).EQ.0).AND.(NSEED.EQ.0).AND.(.NOT.STAY)) & 
                & CALL ACCREJ(NSUCCESS,NFAIL,NSUCCESST,NFAILT)

            TEMP=TEMP*SCALEFAC
            IF (DUMPINT.GT.0) THEN
               IF (MOD(J1,DUMPINT).EQ.0) THEN
                  CALL DUMPSTATE(J1,EBEST,BESTCOORDS,JBEST)
               ENDIF
            ENDIF
           IF (NQ.GT.NSTEPS) GOTO 37
         ENDDO
!  ****************************** END OF LOOP OVER NPAR PARALLEL RUNS *****************************
!
 
         CALL FLUSH(LFH)
      ENDDO
! }}}

37    CONTINUE
         WRITE(LFH,21) NSUCCESST*1.0D0/MAX(1.0D0,1.0D0*(NSUCCESST+NFAILT)),&
     &               STEP,ASTEP,TEMP
21       FORMAT('ACCEPTANCE RATIO FOR RUN=',F12.5,' STEP=',F12.5,' ANGULAR STEP FACTOR=',F12.5,' T=',F12.5)
      RETURN
!OP226>}}} 
      END subroutine

      SUBROUTINE TRANSITION(ENEW,EOLD,ATEST,NP,RANDOM,MCTEMP)
      ! declarations {{{
      ! subroutine 

      DOUBLE PRECISION ENEW, EOLD,  MCTEMP, RANDOM
      LOGICAL ATEST
      INTEGER NP

      ! local

      DOUBLE PRECISION ::   DPRAND, EREF, TEOLD, TENEW, RATIO
      DOUBLE PRECISION TRANS, DISTMIN, DISTMINOLD
      LOGICAL FLAT, evap, evapreject
      INTEGER INDEXOLD, INDEXNEW, J1
      DATA DISTMINOLD /0.0D0/
      COMMON /DMIN/ DISTMIN
      common /ev/ evap, evapreject
      ! }}}
      ! body {{{
         TEOLD=EOLD
         TENEW=ENEW

!  Standard canonical sampling.

         IF (TENEW.LT.TEOLD) THEN
            RANDOM=0.0D0
            ATEST=.TRUE.
         ELSE
            RANDOM=DPRAND()
            IF (DEXP(-(TENEW-TEOLD)/MAX(MCTEMP,1.0D-100)).GT.RANDOM) THEN
               ATEST=.TRUE.
            ELSE
               ATEST=.FALSE.
            ENDIF
         ENDIF

      RETURN 
      ! }}}
      END  SUBROUTINE

      SUBROUTINE ACCREJ(NSUCCESS,NFAIL,NSUCCESST,NFAILT)
      ! declarations {{{
      USE COMMONS

      IMPLICIT NONE

      ! subroutine 
      INTEGER NSUCCESS, NFAIL, NFAILT, NSUCCESST

      ! local 
      INTEGER :: J1, J2
      LOGICAL evap, evapreject
      DOUBLE PRECISION HWMAX,P0,FAC
      common /ev/ evap, evapreject
      ! }}}
      ! body {{{

      P0=1.D0*NSUCCESS/(1.D0*(NSUCCESS+NFAIL))
      
      IF (P0.GT.ACCRAT) THEN
         IF(ARMT) THEN
           FAC=LOG(ARMA*ACCRAT+ARMB)/LOG(ARMA*P0+ARMB)
         ELSE
           FAC=1.05D0
         ENDIF
         IF (FIXBOTH) THEN
         ELSE IF (FIXSTEP) THEN
            IF (.NOT.FIXTEMP) TEMP=TEMP/1.05D0
         ELSE
            STEP=FAC*STEP
            ASTEP=ASTEP*1.05D0
! jwrm2> limit step size for percolation to the cutoff distance for determining connectivity
         ENDIF
      ELSE
         IF(ARMT) THEN
           FAC=LOG(ARMA*ACCRAT+ARMB)/LOG(ARMA*P0+ARMB)
         ELSE
           FAC=1.D0/1.05D0
         ENDIF
         IF (FIXBOTH) THEN
         ELSE IF (FIXSTEP) THEN
            IF (.NOT.FIXTEMP) TEMP=TEMP*1.05D0
         ELSE
            STEP=FAC*STEP
            ASTEP=ASTEP/1.05D0
         ENDIF
      ENDIF
!
! Prevent steps from growing out of bounds. The value of 1000 seems sensible, until
! we do something with such huge dimensions?!
!
      STEP=MIN(STEP,1.0D3)
      OSTEP=MIN(OSTEP,1.0D3)
      ASTEP=MIN(ASTEP,1.0D3)

      WRITE(LFH,'(A,I6,A,F8.4,A,F8.4)') 'Acceptance ratio for previous ',NACCEPT,' steps=',P0,'  FAC=',FAC
      IF (FIXBOTH) THEN
      ELSE IF (FIXSTEP) THEN
         IF(.NOT.FIXTEMP) WRITE(LFH,'(A,F12.4)') 'Temperature is now:',TEMP
      ELSE
         WRITE(LFH,'(A)',ADVANCE='NO') 'Steps are now:'
         WRITE(LFH,'(A,F10.4)',ADVANCE='NO') '  STEP=',STEP    
         IF(ASTEP.GT.0.D0) WRITE(LFH,'(A,F10.4)',ADVANCE='NO')'  ASTEP=',ASTEP 
         IF(.NOT.FIXTEMP) WRITE(LFH,'(A,F10.4)') ' Temperature is now:',TEMP
      ENDIF
!
      NSUCCESST=NSUCCESST+NSUCCESS
      NFAILT=NFAILT+NFAIL
      NSUCCESS=0
      NFAIL=0 
!
      RETURN
      ! }}}
      END SUBROUTINE

      SUBROUTINE MYRESET(NATOMS,NSEED)
      ! {{{
      IMPLICIT NONE

      INTEGER NSEED, J2, NATOMS

      DO J2=1,3*(NATOMS-NSEED)
         COORDS(J2)=COORDSO(J2)
      ENDDO
      DO J2=1,NATOMS
         VAT(J2)=VATO(J2)
      ENDDO

      RETURN
      ! }}}
      END SUBROUTINE
      
      ENDMODULE
