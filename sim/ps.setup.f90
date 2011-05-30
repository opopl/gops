
!> @name SETUP
!> @brief  Read in databases of A and B minima, calculate partition functions, rate constants
!>  and sums of rate constants for all the transition states in the database.
!>
!>  We need pre-existing databases to specify which minima are A and which are B.
!
      SUBROUTINE SETUP
      ! Declarations: Modules and Variables {{{

      USE UTILS
      USE KEY
      USE COMMONS
      USE PORFUNCS
      USE FUNC

      IMPLICIT NONE

      INTEGER J1, J2, STATUS, J3, NDUMMY, NRANDOM, NCOUNT, NMINREMOVE, NTSREMOVE, NMINRETAIN, NTSRETAIN, ISTAT, J4
      DOUBLE PRECISION LOCALPOINTS(3*NATOMS), IXM, IYM, IZM, LOCALPOINTS2(3*NATOMS), DISTANCE, RMAT(3,3), DIST2, DPRAND
      DOUBLE PRECISION PFNORM1, PFNORM2
      DOUBLE PRECISION, ALLOCATABLE :: NEWPFMIN(:)
      INTEGER, ALLOCATABLE :: CANDIDATES(:), MINPREV(:), MINREMOVE(:), TSREMOVE(:), MINRETAIN(:), TSRETAIN(:)
      DOUBLE PRECISION NEWEMIN, NEWIXMIN, NEWIYMIN, NEWIZMIN, NEWFVIBMIN, TSEGUESS, NEWMINCURVE, NEWMINFRQ2,
     &                 TSFVIBGUESS, DUMMY, FRICTIONFAC
      DOUBLE PRECISION :: CUT_UNDERFLOW=-300.0D0
      LOGICAL DEADTS
      INTEGER NEWHORDERMIN!}}}
! subroutine body {{{
      IF (CHARMMT.AND..NOT.MACHINE) CALL READREF('INPUT.CRD')
      
      CALL OPENF(UMIN,"DA",'points.min')
      CALL OPENF(UTS,"DA",'points.ts')

include ps.setup.extractmin.inc.f90
include ps.setup.extractts.inc.f90
include ps.setup.dijinit.inc.f90
include ps.setup.readminA.inc.f90
include ps.setup.readminB.inc.f90
include ps.setup.loadminima.inc.f90

30    IF (YESNO) CLOSE(UMINDATA) ! SAT need to reopen this file

      IF (PRINTT) WRITE(*,'(A,I7,A)') 'setup> parameters read for ',NMIN,' min'
      DO J1=1,NMINA
         IF (LOCATIONA(J1).GT.NMIN) THEN
            PRINT '(3(A,I8))','setup> ERROR - A minimum ',J1,' is number ',LOCATIONA(J1),' but total minima=',NMIN
            STOP
         ENDIF
      ENDDO
      DO J1=1,NMINB
         IF (LOCATIONB(J1).GT.NMIN) THEN
            PRINT '(3(A,I8))','setup> ERROR - B minimum ',J1,' is number ',LOCATIONB(J1),' but total minima=',NMIN
            STOP
         ENDIF
      ENDDO
      IF (CALCORDERT) THEN
         CALL CALCORDER(NATOMS,NMIN,NTS,UMIN,UTS,DEBUG)
         STOP
      ENDIF
      IF (NMIN.GT.0) THEN
         CALL OPENF(UMINDATA,"RW>>","min.data") ! read used in Dijkstra
      ENDIF
      IF ((NMIN.EQ.0).AND.(READMINT.OR.MERGEDBT)) THEN
         CALL OPENF(UMINDATA,">>","min.data") ! for READMIN and MERGEDB startup
      ENDIF
!
!  Check that the necessary coordinates are in fact present. 
!
      !op226{{{
      IF ((.NOT.NOPOINTS).AND.(NATTEMPT.GT.0)) THEN
          IF (AMHT) THEN
            WRITE(*,*) 'setup> AVOIDING MOMENT OF INITERIA CALC FOR AMH'
          ELSE
           DO J1=1,NMIN
            READ(UMIN,REC=J1) (LOCALPOINTS(J2),J2=1,3*NATOMS)
            CALL INERTIAWRAPPER(LOCALPOINTS,NATOMS,angleAxis,IXM,IYM,IZM)
!           IF (PRINTT) WRITE(*,'(2F20.10,I6,3F20.10)') EMIN(J1),FVIBMIN(J1),HORDERMIN(J1),IXM,IYM,IZM
            IF ((ABS(IXM-IXMIN(J1)).GT.IDIFFTOL).OR.&
               (ABS(IYM-IYMIN(J1)).GT.IDIFFTOL).OR.&
               (ABS(IZM-IZMIN(J1)).GT.IDIFFTOL)) THEN
               WRITE(*,'(A,I6)') 'setup> WARNING - principal moments of inertia do not agree with input for minimum ',J1
               WRITE(*,'(A,3F20.10)') 'setup> values from coordinates: ',IXM,IYM,IZM
               WRITE(*,'(A,2F20.10,I6,3F20.10)') 'setup> values from min.data: ', 
     &                     EMIN(J1),FVIBMIN(J1),HORDERMIN(J1),IXMIN(J1),IYMIN(J1),IZMIN(J1)
               PRINT '(A)','LOCALPOINTS:'
               PRINT '(3G20.10)',LOCALPOINTS(1:3*NATOMS)
               IXMIN(J1)=IXM
               IYMIN(J1)=IYM
               IZMIN(J1)=IZM
!              STOP  
            ENDIF
           ENDDO
          ENDIF
         IF (PRINTT) WRITE(*,'(A,I6,A)') 'setup> points for ',NMIN,' minima read from file points.min'
      ENDIF
      !op226 }}}
!

include ps.setup.readmin.inc.f90

      IF (CLOSEFILEST) CLOSE(UNIT=UMINDATA)
!
!  Set all FVIBMIN to ln(2 pi) if NOFRQS is true for consistency.
!  May be useful for running REGROUPFREE on the basis of potential energy only.
!  Won;t work with FREEPAIRS unless the OPTIM runs are also run with NOFRQS.
!
      IF (NOFRQS) FVIBMIN(1:NMIN)=4.675754133D0 ! 2 ln(2pi) +1

! Calculate partition functions for minima
! ====================================
      include ps.setup.calc.pf.inc.f90
! ====================================
! Optional change of reactant minima set via reweighting
! ====================================
      include ps.setup.reweight.inc.f90
! ====================================
! Load transition states
! ====================================
      include ps.setup.read.tsdata.inc.f90 
! ====================================
! Create a ts entry for DUMMYTS runs if there are minima that seem to lack such entries.
! ====================================
      include ps.setup.dummyts.inc.f90

!
!  Set transition state vibrational product to unity for consistency if NOFRQS is true.
!  Won't work with FREEPAIRS unless the OPTIM runs are also run with NOFRQS.
!
      IF (NOFRQS) THEN
         FVIBTS(1:NTS)=1.0D0
         NEGEIG(1:NTS)=-1.0D0
      ENDIF

      IF (CLOSEFILEST) CLOSE(UNIT=UTSDATA)
!
!  Procedure to remove stationary points that are unconnected from A or B sets (or both)
!  according to the prevailing NCONNMIN value.
!
      IF (REMOVEUNCONNECTEDT) THEN
         CALL REMOVE_UNCONNECTED
         STOP
      ENDIF

! ====================================
!  Procedure to remove selected stationary points specified by min.remove and ts.remove.
!  First line of each file gives the numbers of structures to remove.
! ====================================
      include ps.setup.removesp.inc.f90
! ====================================
!  Procedure to retain selected stationary points specified by min.retain.
!  First line of each file gives the numbers of structures to retain.
!  All ts linking minima in the retain list are themselves retained.
! ====================================
      include ps.setup.retainsp.inc.f90
! ====================================
!  Rate constants.
! ====================================
      include ps.setup.rate.const.inc.f90

      IF (MERGEDBT) THEN
         CALL MERGEDB
         STOP
      ENDIF


!  Add transition states and minima from the <PATHNAME> file.
!  Use GETNEWPATH to do the bookkeeping.
! ====================================
        include ps.setup.addpath.inc.f90
! ====================================
        include ps.setup.read.commitdata.inc.f90
! ====================================
!  Read in the pairs of minima previously searched in pairs.data exists.
! ====================================
        include ps.setup.read.pairsdata.inc.f90
! ====================================
!  Read in the minima previously searched in UNTRAP runs.
! ====================================
        include ps.setup.read.mindone.inc.f90
! ====================================
!  Initialise PAIRDIST array for use in making an intial connection.
!  PAIRDIST should contain zero if the two minima are linked by a transition state.
! ====================================
        include ps.setup.init.pairdist.inc.f90
! ====================================
! If USEPAIRST is true then read the sequence of minima from file USEPAIRSFILE
! USEPAIRSFILE must be formatted as a single Epath file
        include ps.setup.usepairs.inc.f90
      
      IF (DOST) CALL DOS
      IF (CVT) CALL CV

      IF ((CONNECTIONS.GT.1).AND.(CHECKCONNECTIONST)) THEN
         WRITE(*,'(A,I6,A)') 'setup> checking for at least ',CONNECTIONS,' connections per minimum'
         DO J1=1,NMIN
            CALL TSSEARCH(J1,0)
         ENDDO
      ENDIF
! }}}
      RETURN 
      END

      DOUBLE PRECISION FUNCTION TSEGUESS(E1,E2,C1,C2,DISTANCE)
      IMPLICIT NONE
      DOUBLE PRECISION E1, E2, DISTANCE, C1, C2, ARGUMENT
!
!  Double fold formulation.!{{{
!      
!     ARGUMENT=c1*c2*distance**2 - 6*(c1 - c2)*(e1 - e2)
!     IF (ARGUMENT.LT.0.0D0) ARGUMENT=0.0D0
!     IF (C1.EQ.C2) THEN
!        TSEGUESS=((c1*distance**2 + 6*e1 - 6*e2)**2/(4.*c1*distance**2) + 6*e2)/6.
!     ELSE
!        TSEGUESS=(c1*c2*(c1 + c2)*distance**2 - 2*c1*c2*distance*Sqrt(ARGUMENT)
!    &             + 6*(c1 - c2)*(-(c2*e1) + c1*e2))/(6.*(c1 - c2)**2)
!     ENDIF
!     IF (TSEGUESS.LT.MAX(E1,E2)) TSEGUESS=MAX(E1,E2)+ABS(E1-E2)!}}}
!
!  Double quadratic formulation.!{{{
!      
!     ARGUMENT=c1*c2*distance**2 - 2*(c1 - c2)*(e1 - e2)
!     IF (ARGUMENT.LT.0.0D0) ARGUMENT=0.0D0
!     IF (C1.EQ.C2) THEN
!        TSEGUESS=((c1*distance**2)/4. + e1 + (e1 - e2)**2/(c1*distance**2) + e2)/2.
!     ELSE
!        TSEGUESS=(c1*c2*(c1 + c2)*distance**2 - 2*c1*c2*distance*
!    &     Sqrt(ARGUMENT) + 2*(c1 - c2)*(-(c2*e1) + c1*e2))/(2.*(c1 - c2)**2)
!     ENDIF!}}}

      TSEGUESS=MAX(E1,E2)+DISTANCE
      
      END FUNCTION TSEGUESS

      DOUBLE PRECISION FUNCTION TSFVIBGUESS(E1,E2,FVIB1,FVIB2,MINF1,MINF2,NATOMS)
      IMPLICIT NONE
      DOUBLE PRECISION E1, E2, FVIB1,FVIB2, MINF1, MINF2
      INTEGER NATOMS
!
!  The conversion factor for CHARMM and AMBER is included in the MINFRQ2 values read from min.data.info
!  The MINFRQ2 values are read as the ln from min.data.info
!
!     IF (E1.GT.E2) THEN
!        TSFVIBGUESS=FVIB1-MINF1
!     ELSE
!        TSFVIBGUESS=FVIB2-MINF2
!     ENDIF
      IF (E1.GT.E2) THEN
         TSFVIBGUESS=FVIB1*(3*NATOMS-7)/(3*NATOMS-6)
      ELSE
         TSFVIBGUESS=FVIB2*(3*NATOMS-7)/(3*NATOMS-6)
      ENDIF

      END FUNCTION TSFVIBGUESS
