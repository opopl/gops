
!>  Read in databases of A and B minima, calculate partition functions, rate constants
!>  and sums of rate constants for all the transition states in the database.
!
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

!  In this case we try to connect two minima, whose coordinates are in
!  files odata.start and odata.finish for a DIJINITSTARTT run. We set up min.data as 
!  well as min.A and min.B, with entries 1 1 and 1 2.
!
!  DIJINITCONTT specifies a continuation of an initial path run.
!
      !op226 IF (DIJINITT.OR.DIJINITFLYT) THEN {{{
      IF (DIJINITT.OR.DIJINITFLYT) THEN
         INQUIRE(FILE='min.data',EXIST=YESNO)
         IF (DIJINITCONTT) THEN
            IF (.NOT.YESNO) THEN
               PRINT '(A)','setup> ERROR - min.data must exist for a DIJINIT continuation run'
               STOP
            ENDIF
         ELSEIF (DIJINITSTARTT) THEN
            IF (YESNO) CALL MYSYSTEM(STATUS,DEBUG,'cp min.data min.data.save')
            PRINT '(A)','setup> Creating new min.data file'
            OPEN(UNIT=UMINDATA,FILE='min.data',STATUS='UNKNOWN')
            INQUIRE(FILE='pairs.data',EXIST=YESNO)
            IF (YESNO) THEN
               PRINT '(A)','setup> Moving old pairs.data file'
               CALL MYSYSTEM(STATUS,DEBUG,'mv pairs.data pairs.data.old')
            ENDIF
            INQUIRE(FILE='commit.data',EXIST=YESNO)
            IF (YESNO) THEN
               PRINT '(A)','setup> Moving old commit.data file'
               CALL MYSYSTEM(STATUS,DEBUG,'mv commit.data commit.data.old')
            ENDIF
            INQUIRE(FILE='min.data.info.1',EXIST=YESNO)
            IF (YESNO) CALL MYSYSTEM(STATUS,DEBUG,'rm min.data.info.1')
            CALL MYSYSTEM(STATUS,DEBUG,'cp odata.start odata.1')
            IF (CHARMMT) CALL MYSYSTEM(STATUS,DEBUG,'cp start.crd input.crd.1')
            IF (AMHT) CALL MYSYSTEM(STATUS,DEBUG,'cp AMHstart start.1')
            CALL MYSYSTEM(STATUS,DEBUG,TRIM(ADJUSTL(EXEC))//' 1 > output.start')
            INQUIRE(FILE='min.data.info.1',EXIST=YESNO)
            IF (YESNO) THEN
               OPEN(UNIT=1,FILE='min.data.info.1',STATUS='OLD')
               IF (DUMMYTST.AND.LOWESTFRQT) THEN
                  READ(1,*) EMIN(1),FVIBMIN(1),HORDERMIN(1),IXMIN(1),IYMIN(1),IZMIN(1),MINCURVE(1),MINFRQ2(1)
                  WRITE(UMINDATA,'(2F20.10,I6,5F20.10)') EMIN(1), FVIBMIN(1), HORDERMIN(1),IXMIN(1),
     &                                                      IYMIN(1),IZMIN(1),MINCURVE(1),MINFRQ2(1)
               ELSE
                  READ(1,*) EMIN(1),FVIBMIN(1),HORDERMIN(1),IXMIN(1),IYMIN(1),IZMIN(1)
                  WRITE(UMINDATA,'(2F20.10,I6,3F20.10)') EMIN(1), FVIBMIN(1), HORDERMIN(1),IXMIN(1), IYMIN(1), IZMIN(1)
               ENDIF
               READ(1,*) (LOCALPOINTS(J2),J2=1,3*NATOMS)
               WRITE(UMIN,REC=1) (LOCALPOINTS(J2),J2=1,3*NATOMS)
               CLOSE(1)
            ELSE
               PRINT *, 'setup> ERROR - no min.data.info.1 found - check OPTIM output in output.start'
               STOP        
            ENDIF
            INQUIRE(FILE='min.data.info.2',EXIST=YESNO)
            IF (YESNO) CALL MYSYSTEM(STATUS,DEBUG,'rm min.data.info.2')
            CALL MYSYSTEM(STATUS,DEBUG,'cp odata.finish odata.2')
            IF (CHARMMT) CALL MYSYSTEM(STATUS,DEBUG,'cp finish.crd input.crd.2')
            IF (AMHT) CALL MYSYSTEM(STATUS,DEBUG,'cp AMHfinish start.2')
            CALL MYSYSTEM(STATUS,DEBUG,TRIM(ADJUSTL(EXEC))//' 2 > output.finish')
            INQUIRE(FILE='min.data.info.2',EXIST=YESNO)
            IF (YESNO) THEN
               OPEN(UNIT=1,FILE='min.data.info.2',STATUS='OLD')
               IF (DUMMYTST.AND.LOWESTFRQT) THEN
                  READ(1,*) EMIN(2),FVIBMIN(2),HORDERMIN(2),IXMIN(2),IYMIN(2),IZMIN(2),MINCURVE(2),MINFRQ2(2)
                  WRITE(UMINDATA,'(2F20.10,I6,5F20.10)') EMIN(2), FVIBMIN(2), HORDERMIN(2),IXMIN(2), 
     &                                                   IYMIN(2),IZMIN(2),MINCURVE(2),MINFRQ2(2)
               ELSE
                  READ(1,*) EMIN(2),FVIBMIN(2),HORDERMIN(2),IXMIN(2),IYMIN(2),IZMIN(2)
                  WRITE(UMINDATA,'(2F20.10,I6,3F20.10)') EMIN(2), FVIBMIN(2), HORDERMIN(2),IXMIN(2), IYMIN(2), IZMIN(2)
               ENDIF
               CLOSE(UMINDATA)
               READ(1,*) (LOCALPOINTS(J2),J2=1,3*NATOMS)
               WRITE(UMIN,REC=2) (LOCALPOINTS(J2),J2=1,3*NATOMS)
               CLOSE(1)
            ELSE
               PRINT *, 'setup> ERROR - no min.data.info.2 found - check OPTIM output in output.finish'
               STOP
            ENDIF
            INQUIRE(FILE='min.A',EXIST=YESNO)
            IF (DIJINITSTARTT) THEN
               IF (YESNO) CALL MYSYSTEM(STATUS,DEBUG,'cp min.A min.A.save')
               PRINT '(A)','setup> Creating new min.A file'
            ELSEIF (DIJINITCONTT) THEN
               IF (.NOT.YESNO) THEN
                  PRINT '(A)','setup> ERROR - min.A must exist for a DIJINIT continuation run'
                  STOP
               ENDIF
            ENDIF
            OPEN(1,FILE='min.A',STATUS='UNKNOWN')
            WRITE(1,'(I6)') 1
            WRITE(1,'(I6)') 1
            CLOSE(1)
            INQUIRE(FILE='min.B',EXIST=YESNO)
            IF (DIJINITSTARTT) THEN
               IF (YESNO) CALL MYSYSTEM(STATUS,DEBUG,'cp min.B min.B.save')
               PRINT '(A)','setup> Creating new min.B file'
            ELSEIF (DIJINITCONTT) THEN
               IF (.NOT.YESNO) THEN
                  PRINT '(A)','setup> ERROR - min.B must exist for a DIJINIT continuation run'
                  STOP
               ENDIF
            ENDIF
            OPEN(1,FILE='min.B',STATUS='UNKNOWN')
            WRITE(1,'(I6)') 1
            WRITE(1,'(I6)') 2
            CLOSE(1)
            PRINT '(A)','setup> initial OPTIM jobs run for odata.start and odata.finish'

            IF (DUMMYTST) THEN
               READ(UMIN,REC=1) (LOCALPOINTS(J2),J2=1,3*NATOMS)
               READ(UMIN,REC=2) (LOCALPOINTS2(J3),J3=1,3*NATOMS)
               CALL MINPERMDIST(LOCALPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY,
     &                          RMAT,.FALSE.)
               IF (INTERPCOSTFUNCTION) CALL MINPERMDIST(LOCALPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD, 
     &                                                  DISTANCE,DIST2,RIGIDBODY,RMAT,INTERPCOSTFUNCTION)
               MINDISTMIN(1)=DISTANCE
               MINDISTMIN(2)=DISTANCE
!
!  Must create an entry in ts.data in this case.
!  ETS,FVIBTS,HORDERTS,PLUS,MINUS,IXTS,IYTS,IZTS
!
               INQUIRE(FILE='ts.data',EXIST=YESNO)
               IF (YESNO) THEN
                  PRINT '(A)','setup> - file ts.data already exists. Copying to ts.data.save'
                  CALL MYSYSTEM(STATUS,DEBUG,'mv ts.data ts.data.save')
               ENDIF
               OPEN(UNIT=UTSDATA,FILE='ts.data',STATUS='NEW') 

               IF (IMFRQT) THEN
                  PRINT '()',"setup> ERROR: can''t guess negative eigenvalue - don''t use DUMMYTS and IMFRQ"
                  STOP
               ENDIF
               WRITE(UTSDATA,'(2F20.10,3I10,3F20.10)') TSEGUESS(EMIN(1),EMIN(2),MINCURVE(1),MINCURVE(2),DISTANCE), 
     &               TSFVIBGUESS(EMIN(1),EMIN(2),FVIBMIN(1),FVIBMIN(2),MINFRQ2(1),MINFRQ2(2),NATOMS),1,1,2,1.0D0,1.0D0,1.0D0

               CLOSE(UTSDATA)
            ENDIF
         ENDIF
      ENDIF
      !op226 }}}

      !op226 IF (STARTFROMPATH) THEN {{{

      IF (STARTFROMPATH) THEN
!
!  Get all the necessary information about the A and B minima from the <PATHNAME> file.
!  The assumption is that we have just one A minimum and one B in this case.
!  Use GETNEWPATH or GETALLPATHS to do the bookkeeping.
!  Detect presence of existing min.data file and STOP if found to prevent overwriting.
!
         INQUIRE(FILE='min.data',EXIST=YESNO)
         IF (YESNO) THEN
            PRINT '(A)','ERROR - file min.data already exists. Will not overwrite.'
            STOP
         ENDIF
         CALL MYSYSTEM(STATUS,DEBUG,'cp ' // TRIM(ADJUSTL(PATHNAME)) // ' path.info')
         OPEN(UNIT=UMINDATA,FILE='min.data',STATUS='NEW')
         INQUIRE(FILE='ts.data',EXIST=YESNO)
         IF (YESNO) THEN
            PRINT '(A)','ERROR - file ts.data already exists. Will not overwrite.'
            STOP
         ENDIF
         OPEN(UNIT=UTSDATA,FILE='ts.data',STATUS='NEW')
         IF (CLOSEFILEST) CLOSE(UNIT=UTSDATA)
         IF (CLOSEFILEST) CLOSE(UNIT=UMINDATA)
         IF (STARTTRIPLES) THEN
            CALL GETALLPATHS
         ELSE
            CALL GETNEWPATH(0,0)
         ENDIF
         INQUIRE(FILE='min.A',EXIST=YESNO)
         IF (YESNO) THEN
            PRINT '(A)','ERROR - file min.A already exists. Will not overwrite.'
            STOP
         ENDIF
         OPEN(1,FILE='min.A',STATUS='NEW') 
         WRITE(1,'(I6)') 1
         WRITE(1,'(I6)') STARTMINA
         CLOSE(1)
         INQUIRE(FILE='min.B',EXIST=YESNO)
         IF (YESNO) THEN
            PRINT '(A)','ERROR - file min.B already exists. Will not overwrite.'
            STOP
         ENDIF
         OPEN(1,FILE='min.B',STATUS='NEW') 
         WRITE(1,'(I6)') 1
         WRITE(1,'(I6)') STARTMINB
         CLOSE(1)
         CLOSE(UMINDATA)
         CLOSE(UTSDATA)
         PRINT '(A,I5,A,I5,A)','setup> The unique A and B minima are ',STARTMINA,' and ',STARTMINB,' respectively'
      ENDIF
      !op226 }}}
   
      !op226 open min.A; read in locations of minima {{{
      
      INQUIRE(FILE='min.A',EXIST=YESNO)
      IF (YESNO) THEN
         OPEN(1,FILE='min.A',STATUS='OLD')
         READ(1,*) NMINA
         ALLOCATE(LOCATIONA(NMINA)) ! not deallocated
         READ(1,*) LOCATIONA(1:NMINA)
         CLOSE(1)
      ELSEIF (.NOT.(READMINT.OR.MERGEDBT)) THEN
         PRINT '(A)','setup> ERROR - no min.A file'
         STOP
      ELSE
         NMINA=0
      ENDIF
      IF (DEBUG) THEN
         PRINT '(A,I6,A)','setup> there are ',NMINA,' A minima at locations:'
         PRINT '(10I6)',LOCATIONA(1:NMINA)
      ENDIF
      IF (PRINTT) WRITE(*,'(A,I6,A)') 'setup> locations read for ',NMINA,' min of type A'

      !op226 }}}

      !op226 open min.B; read in locations of minima {{{

      INQUIRE(FILE='min.B',EXIST=YESNO)
      IF (YESNO) THEN
         OPEN(1,FILE='min.B',STATUS='OLD')
         READ(1,*) NMINB
         ALLOCATE(LOCATIONB(NMINB)) ! not deallocated
         READ(1,*) LOCATIONB(1:NMINB)
         CLOSE(1)
      ELSEIF (.NOT.(READMINT.OR.MERGEDBT)) THEN
         PRINT '(A)','setup> ERROR - no min.B file'
         STOP
      ELSE
         NMINB=0
      ENDIF
      IF (DEBUG) THEN
         PRINT '(A,I6,A)','setup> there are ',NMINB,' B minima at locations:'
         PRINT '(10I6)',LOCATIONB(1:NMINB)
      ENDIF
      IF (PRINTT) WRITE(*,'(A,I6,A)') 'setup> locations read for ',NMINB,' min of type B'

      !op226 }}}
!
!  Load the minima.
!
      !op226 open min.data; load the minima {{{

      INQUIRE(FILE='min.data',EXIST=YESNO)
      IF (YESNO) THEN
         OPEN(UNIT=UMINDATA,FILE='min.data',STATUS='OLD')
         NMIN=0
         J1=0
         DO 
            J1=J1+1
            IF (J1.GT.MAXMIN) CALL MINDOUBLE
            IF (DUMMYTST.AND.LOWESTFRQT) THEN
               READ(UMINDATA,*,END=30) EMIN(J1),FVIBMIN(J1),HORDERMIN(J1),IXMIN(J1),IYMIN(J1),IZMIN(J1),MINCURVE(J1),MINFRQ2(J1)
            ELSE
               READ(UMINDATA,*,END=30) EMIN(J1),FVIBMIN(J1),HORDERMIN(J1),IXMIN(J1),IYMIN(J1),IZMIN(J1)
            ENDIF
            NMIN=NMIN+1
         ENDDO
      ELSEIF (.NOT.(READMINT.OR.MERGEDBT)) THEN
         PRINT '(A)','setup> ERROR - no min.data file'
         STOP
      ELSE
         NMIN=0
      ENDIF

      !op226 }}}
      
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
         CALL OPENF(UMINDATA,">>","min.data")
         !OPEN(UNIT=UMINDATA,FILE='min.data',STATUS='OLD',POSITION="APPEND",ACTION="READWRITE") ! read used in Dijkstra
      ENDIF
      IF ((NMIN.EQ.0).AND.(READMINT.OR.MERGEDBT)) THEN
         OPEN(UNIT=UMINDATA,FILE='min.data',STATUS='UNKNOWN',POSITION="APPEND") ! for READMIN and MERGEDB startup
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
            IF ((ABS(IXM-IXMIN(J1)).GT.IDIFFTOL).OR.
     1          (ABS(IYM-IYMIN(J1)).GT.IDIFFTOL).OR.
     1          (ABS(IZM-IZMIN(J1)).GT.IDIFFTOL)) THEN
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
! Read data for minima from min.data.info style file. Multiple minima are allowed.
!
!op226 IF (READMINT) THEN {{{
      IF (READMINT) THEN
         INQUIRE(FILE=TRIM(ADJUSTL(MINNAME)),EXIST=YESNO)
         IF (YESNO) THEN
            OPEN(UNIT=1,FILE=TRIM(ADJUSTL(MINNAME)),STATUS='OLD')
            DO 
               IF (DUMMYTST.AND.LOWESTFRQT) THEN
                  READ(1,*,END=130) NEWEMIN,NEWFVIBMIN,NEWHORDERMIN,NEWIXMIN,NEWIYMIN,NEWIZMIN,NEWMINCURVE,NEWMINFRQ2
               ELSE
                  READ(1,*,END=130) NEWEMIN,NEWFVIBMIN,NEWHORDERMIN,NEWIXMIN,NEWIYMIN,NEWIZMIN
               ENDIF
               READ(1,*) (LOCALPOINTS(J2),J2=1,3*NATOMS)
!
! Must check it is not an old minimum!
!
               DO J2=1,NMIN
                  DISTANCE=1.0D100
                  IF (ABS(NEWEMIN-EMIN(J2)).LT.EDIFFTOL) THEN
                     READ(UMIN,REC=J2) (LOCALPOINTS2(J3),J3=1,3*NATOMS)
                     CALL MINPERMDIST(LOCALPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE, 
     &                                DIST2,RIGIDBODY,RMAT,.FALSE.)
                  ENDIF
      
                  IF ((ABS(NEWEMIN-EMIN(J2)).LT.EDIFFTOL).AND.(DISTANCE.LT.GEOMDIFFTOL)) THEN
                     PRINT '(A,I6)','setup> minimum is database minimum ',J2
                     IF (ABS(NEWFVIBMIN-FVIBMIN(J2))/FVIBMIN(J2).GT.1.0D-4) THEN
                        WRITE(*,'(A,F15.5,A,F15.5)') 'setup> WARNING, NEWFVIBMIN=',NEWFVIBMIN,' should be ',FVIBMIN(J2)
                     ENDIF
                     IF (NEWHORDERMIN.NE.HORDERMIN(J2)) THEN
                        WRITE(*,'(A,I6,A,I6)') 'setup> ERROR, NEWHORDERMIN=',NEWHORDERMIN,' should be ',HORDERMIN(J2)
                        NEWHORDERMIN=MAX(NEWHORDERMIN,HORDERMIN(J2))
                        WRITE(*,'(A,I6)') 'setup> using maximum value: ',NEWHORDERMIN
                     ENDIF
                     GOTO 140
                  ENDIF
               ENDDO
               NMIN=NMIN+1
               IF (NMIN.GT.MAXMIN) CALL MINDOUBLE
               EMIN(NMIN)=NEWEMIN
               FVIBMIN(NMIN)=NEWFVIBMIN
               HORDERMIN(NMIN)=NEWHORDERMIN
               IXMIN(NMIN)=NEWIXMIN
               IYMIN(NMIN)=NEWIYMIN
               IZMIN(NMIN)=NEWIZMIN
               IF (DUMMYTST.AND.LOWESTFRQT) MINCURVE(NMIN)=NEWMINCURVE
               IF (DUMMYTST.AND.LOWESTFRQT) MINFRQ2(NMIN)=NEWMINFRQ2
               WRITE(*,'(A,I6,A)') 'setup> new minimum ',NMIN,
     &                ' writing parameters to file min.data and points to points.min'
               IF (DUMMYTST.AND.LOWESTFRQT) THEN
                  WRITE(UMINDATA,'(2F20.10,I6,5F20.10)') EMIN(NMIN), FVIBMIN(NMIN), HORDERMIN(NMIN), 
     &                                             IXMIN(NMIN), IYMIN(NMIN), IZMIN(NMIN), MINCURVE(NMIN),MINFRQ2(NMIN)
               ELSE
                  WRITE(UMINDATA,'(2F20.10,I6,3F20.10)') EMIN(NMIN), FVIBMIN(NMIN), HORDERMIN(NMIN), 
     &                                             IXMIN(NMIN), IYMIN(NMIN), IZMIN(NMIN)
               ENDIF
               CALL FLUSH(UMINDATA,ISTAT)
               WRITE(UMIN,REC=NMIN) (LOCALPOINTS(J2),J2=1,3*NATOMS)
140            CONTINUE
            ENDDO
130         CLOSE(1)
         ELSE
            PRINT '(A)','setup> ERROR - no file ',TRIM(ADJUSTL(MINNAME))
         ENDIF
         STOP
      ENDIF
      !op226 }}}

      IF (CLOSEFILEST) CLOSE(UNIT=UMINDATA)
!
!  Set all FVIBMIN to ln(2 pi) if NOFRQS is true for consistency.
!  May be useful for running REGROUPFREE on the basis of potential energy only.
!  Won;t work with FREEPAIRS unless the OPTIM runs are also run with NOFRQS.
!
      IF (NOFRQS) FVIBMIN(1:NMIN)=4.675754133D0 ! 2 ln(2pi) +1
!
!  Calculate partition functions for minima. Note that the total partition function
!  is not needed, just the totals for A and B. Since A and B sets are fixed here
!  we don;t need to change the totals.
!
!     PFMEAN=0.0D0
      PFMEAN=-HUGE(1.0D0)
!     PFNORM1=0.0D0 ! use this to calculate ratios without the pe factor
!     PFNORM2=0.0D0 ! use this to calculate ratios with the pe factor

      !op226 {{{

      IF (ENSEMBLE.EQ.'T') THEN
         IF (TEMPERATURE.LE.0.0D0) THEN
            PRINT '(A)','setup> ERROR - TEMPERATURE=',TEMPERATURE
            STOP
         ENDIF
         DO J1 = 1,NMIN
            PFMIN(J1) = -EMIN(J1)/TEMPERATURE - FVIBMIN(J1)/2.0D0 - LOG(1.0D0*HORDERMIN(J1))
!           PFMEAN=PFMEAN+PFMIN(J1)
!           PFNORM1=PFNORM1+EXP(- FVIBMIN(J1)/2.0D0 - LOG(1.0D0*HORDERMIN(J1)))
!           PFNORM2=PFNORM2+EXP(PFMIN(J1))
            IF (PFMIN(J1).GT.PFMEAN) PFMEAN=PFMIN(J1)
         ENDDO
      ELSEIF (ENSEMBLE.EQ.'E') THEN
         DO J1 = 1,NMIN
            IF (TOTALE.GT.EMIN(J1)) THEN
               PFMIN(J1) = (KAPPA-1)*LOG(TOTALE-EMIN(J1)) - FVIBMIN(J1)/2.0D0 - LOG(1.0D0*HORDERMIN(J1))
!              PFMEAN=PFMEAN+PFMIN(J1)
!              PFNORM1=PFNORM1+EXP(- FVIBMIN(J1)/2.0D0 - LOG(1.0D0*HORDERMIN(J1)))
!              PFNORM2=PFNORM2+EXP(PFMIN(J1))
               IF (PFMIN(J1).GT.PFMEAN) PFMEAN=PFMIN(J1)
            ELSE
               PFMIN(J1) = -1.0D250
            ENDIF
         ENDDO
      ELSE
         PRINT*,'ERROR, ENSEMBLE must be set to T or E'
         STOP
      ENDIF
      !op226 }}}
!     PFMEAN=PFMEAN/NMIN ! DJW
      IF (DEBUG) THEN
         WRITE(*,'(A,G20.10)') 'setup> mean ln Z=',PFMEAN
!        WRITE(*,'(A)') '     energy        pg order     high T/E prob       Peq'
!        DO J1=1,NMIN
!           WRITE(*,'(F20.10,I6,2G20.10)') EMIN(J1),HORDERMIN(J1), 
!    &                    EXP(-FVIBMIN(J1)/2.0D0-LOG(1.0D0*HORDERMIN(J1))-LOG(PFNORM1)), 
!    &                    EXP(PFMIN(J1)-LOG(PFNORM2))
!        ENDDO
      ENDIF
      DO J1=1,NMIN
         PFMIN(J1) = PFMIN(J1) - PFMEAN
      ENDDO

      PFTOTALB=0.0D0
      DO J1=1,NMINB
         PFTOTALB=PFTOTALB+EXP(PFMIN(LOCATIONB(J1))-PFMIN(LOCATIONB(1)))
      ENDDO
      IF (NMINB.GT.0.0D0) PFTOTALB=LOG(PFTOTALB)+PFMIN(LOCATIONB(1))

      PFTOTALA=0.0D0
      DO J1=1,NMINA
         PFTOTALA=PFTOTALA+EXP(PFMIN(LOCATIONA(J1))-PFMIN(LOCATIONA(1)))
      ENDDO
      IF (NMINA.GT.0.0D0) PFTOTALA=LOG(PFTOTALA)+PFMIN(LOCATIONA(1))
!
!     Optional change of reactant minima set via reweighting.
!
      IF (REWEIGHTT) THEN
         ALLOCATE(CANDIDATES(NMIN))
         IF (DIRECTION.EQ.'AB') THEN
            ALLOCATE(NEWPFMIN(NMINB))
            NEWPFMIN(1:NMINB)=0.0D0
!
!  Select NRWREACTANT minima from the B set according to the required weights in RWPROB
!
            PFTOTALB=0.0D0
            DO J1=1,NRWBINS
               IF (NINT(NRWREACTANT*RWPROB(J1)).EQ.0) CYCLE
               NCOUNT=0
               DO J2=1,NMINB ! identify minima in the required energy range
                  IF ((EMIN(LOCATIONB(J2)).GE.RWEMIN+(J1-1)*RWBINWIDTH).AND.
     &                (EMIN(LOCATIONB(J2)).LE.RWEMIN+J1*RWBINWIDTH)) THEN
                      NCOUNT=NCOUNT+1
                      CANDIDATES(NCOUNT)=J2
                  ENDIF
               ENDDO
               PRINT '(3(A,I8),A,G20.10)','setup> number of B minima in energy bin ',J1,' is ',NCOUNT,' number needed=',
     &                           NINT(NRWREACTANT*RWPROB(J1)),' probability=',RWPROB(J1)
               IF (NCOUNT.EQ.0) STOP
               DO J2=1,NINT(NRWREACTANT*RWPROB(J1))
                  NRANDOM=NCOUNT*DPRAND()+1
                  PRINT '(3(A,I8))','setup> selecting B minimum number ',CANDIDATES(NRANDOM),
     &                              ' location ',LOCATIONB(CANDIDATES(NRANDOM)),' for the product set'
                  NEWPFMIN(CANDIDATES(NRANDOM))=NEWPFMIN(CANDIDATES(NRANDOM))+1.0D0
               ENDDO
               PFTOTALB=PFTOTALB+NINT(NRWREACTANT*RWPROB(J1))
            ENDDO
            PFTOTALB=LOG(PFTOTALB) ! partition functions are stored as log's
            NCOUNT=0
            DO J1=1,NMINB
               IF (NEWPFMIN(J1).GT.0.0D0) THEN
                  NCOUNT=NCOUNT+1
                  LOCATIONB(NCOUNT)=LOCATIONB(J1)
                  PFMIN(LOCATIONB(NCOUNT))=LOG(NEWPFMIN(J1))
                  PRINT '(A,I8,A,G20.10)','setup> relative weight for reactant minimum ',LOCATIONB(NCOUNT),' is ',
     &                        EXP(PFMIN(LOCATIONB(NCOUNT))-PFTOTALB)
               ENDIF
            ENDDO
            NMINB=NCOUNT
            PRINT '(A,I8,A)','setup> there are now ',NMINB,' minima of type B'
         ELSE
            ALLOCATE(NEWPFMIN(NMINA))
            NEWPFMIN(1:NMINA)=0.0D0
!
!  Select NRWREACTANT minima from the A set according to the required weights in RWPROB
!
            PFTOTALA=0.0D0
            DO J1=1,NRWBINS
               IF (NINT(NRWREACTANT*RWPROB(J1)).EQ.0) CYCLE
               NCOUNT=0
               DO J2=1,NMINA ! identify minima in the required energy range
                  IF ((EMIN(LOCATIONA(J2)).GE.RWEMIN+(J1-1)*RWBINWIDTH).AND.
     &                (EMIN(LOCATIONA(J2)).LE.RWEMIN+J1*RWBINWIDTH)) THEN
                      NCOUNT=NCOUNT+1
                      CANDIDATES(NCOUNT)=J2
                  ENDIF
               ENDDO
               PRINT '(3(A,I8),A,G20.10)','setup> number of A minima in energy bin ',J1,' is ',NCOUNT,' number needed=',
     &                           NINT(NRWREACTANT*RWPROB(J1)),' probability=',RWPROB(J1)
               IF (NCOUNT.EQ.0) STOP
               DO J2=1,NINT(NRWREACTANT*RWPROB(J1))
                  NRANDOM=NCOUNT*DPRAND()+1
                  PRINT '(3(A,I8))','setup> selecting A minimum number ',CANDIDATES(NRANDOM),
     &                              ' location ',LOCATIONA(CANDIDATES(NRANDOM)),' for the product set'
                  NEWPFMIN(CANDIDATES(NRANDOM))=NEWPFMIN(CANDIDATES(NRANDOM))+1.0D0
               ENDDO
               PFTOTALA=PFTOTALA+NINT(NRWREACTANT*RWPROB(J1))
            ENDDO
            PFTOTALA=LOG(PFTOTALA) ! partition functions are stored as log's
            NCOUNT=0
            DO J1=1,NMINA
               IF (NEWPFMIN(J1).GT.0.0D0) THEN
                  NCOUNT=NCOUNT+1
                  LOCATIONA(NCOUNT)=LOCATIONA(J1)
                  PFMIN(LOCATIONA(NCOUNT))=LOG(NEWPFMIN(J1))
                  PRINT '(A,I8,A,G20.10)','setup> relative weight for reactant minimum ',LOCATIONA(NCOUNT),' is ',
     &                        EXP(PFMIN(LOCATIONA(NCOUNT))-PFTOTALA)
               ENDIF
            ENDDO
            NMINA=NCOUNT
            PRINT '(A,I8,A)','setup> there are now ',NMINA,' minima of type A'
         ENDIF
         DEALLOCATE(NEWPFMIN,CANDIDATES)
      ENDIF
!
!  Load transition states.
!
      DO J1=1,NMIN
         TOPPOINTER(J1)=-1
      ENDDO
      INQUIRE(FILE='ts.data',EXIST=YESNO)
      IF (YESNO.AND.DIJINITSTARTT) THEN
         IF (.NOT.DUMMYTST) THEN
            CALL MYSYSTEM(STATUS,DEBUG,'mv ts.data ts.data.save')
            PRINT '(A)','setup> Removing old ts.data file'
            YESNO=.FALSE.
         ENDIF
      ENDIF
      
      IF (YESNO) THEN
         OPEN(UNIT=UTSDATA,FILE='ts.data',STATUS='OLD')
         J1=0
         DO 
            J1=J1+1
            IF (J1.GT.MAXTS) CALL TSDOUBLE
            IF (IMFRQT) THEN
               READ(UTSDATA,*,END=40) ETS(J1),FVIBTS(J1),HORDERTS(J1),PLUS(J1),MINUS(J1),IXTS(J1),IYTS(J1),IZTS(J1),NEGEIG(J1)
            ELSE
               READ(UTSDATA,*,END=40) ETS(J1),FVIBTS(J1),HORDERTS(J1),PLUS(J1),MINUS(J1),IXTS(J1),IYTS(J1),IZTS(J1)
            ENDIF
            IF ((PLUS(J1).GT.NMIN).OR.(MINUS(J1).GT.NMIN)) THEN
               PRINT '(A,I6,A)','setup> ERROR - minima specified for ts ',J1,' lie beyond those specified in min.data'
               PRINT '(A,2I6)','setup> plus and minus minima are ',PLUS(J1),MINUS(J1)
               STOP
            ENDIF

            IF (DUMMYTST.AND.(.NOT.NOPOINTS).AND.(NATTEMPT.GT.0)) THEN
               READ(UMIN,REC=PLUS(J1)) (LOCALPOINTS(J2),J2=1,3*NATOMS)
               READ(UMIN,REC=MINUS(J1)) (LOCALPOINTS2(J3),J3=1,3*NATOMS)
                  CALL MINPERMDIST(LOCALPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY,
     &                             RMAT,.FALSE.)
                  IF (INTERPCOSTFUNCTION) CALL MINPERMDIST(LOCALPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD, 
     &                                                     DISTANCE,DIST2,RIGIDBODY,RMAT,INTERPCOSTFUNCTION)
               IF (DISTANCE.LT.MINDISTMIN(PLUS(J1))) MINDISTMIN(PLUS(J1))=DISTANCE
               IF (DISTANCE.LT.MINDISTMIN(MINUS(J1))) MINDISTMIN(MINUS(J1))=DISTANCE
            ENDIF
         ENDDO
40       CLOSE(UTSDATA) ! SAT need to reopen this file
         OPEN(UNIT=UTSDATA,FILE='ts.data',STATUS='OLD',POSITION="APPEND",ACTION="READWRITE",FORM="FORMATTED") ! read used in Dijkstra
         NTS=J1-1
         IF (PRINTT) WRITE(*,'(A,I7,A)') 'setup> parameters read for ',NTS,' ts'

         IF (DIJKSTRAT .OR. KSHORTESTPATHST) THEN
            INQUIRE(FILE='ts.attempts',EXIST=YESNO)
            TSATTEMPT(1:NTS)=0
            IF (YESNO) THEN
               PRINT '(A)','setup> Reading the number of searches for existing transition states from ts.attempts'
               OPEN(UNIT=1,FILE='ts.attempts',STATUS='UNKNOWN')
               J2=0
               DO J1=1,NTS
                  READ(1,'(I8)',END=51) TSATTEMPT(J1)
                  J2=J2+1
               ENDDO
51             CLOSE(1)
               IF (J2.LT.NTS) PRINT '(A)','setup> WARNING - end of file ts.attempts, remaining attempts set to zero'
            ENDIF
         ENDIF

         DO J1=1,NTS
            POINTERP(J1)=-1
            POINTERM(J1)=-1
         ENDDO
         DO J1=NTS,1,-1
            IF (J1.GT.TOPPOINTER(PLUS(J1)))  TOPPOINTER(PLUS(J1))=J1
            IF (J1.GT.TOPPOINTER(MINUS(J1))) TOPPOINTER(MINUS(J1))=J1
            DO J2=J1-1,1,-1
               IF (PLUS(J2).EQ.PLUS(J1)) THEN
                  POINTERP(J1)=J2
                  GOTO 41
               ELSE IF (MINUS(J2).EQ.PLUS(J1)) THEN
                  POINTERP(J1)=J2
                  GOTO 41
               ENDIF
            ENDDO
41          CONTINUE
            DO J2=J1-1,1,-1
               IF (PLUS(J2).EQ.MINUS(J1)) THEN
                  POINTERM(J1)=J2
                  GOTO 42
               ELSE IF (MINUS(J2).EQ.MINUS(J1)) THEN
                  POINTERM(J1)=J2
                  GOTO 42
               ENDIF
            ENDDO
42          CONTINUE
         ENDDO
!        IF (DEBUG) THEN
!           DO J1=1,NMIN
!              WRITE(*,'(A,I6,A,I6)') 'setup> for minimum ',J1,' last occurrence is for ts number ',TOPPOINTER(J1)
!           ENDDO
!           DO J1=1,NTS
!              WRITE(*,'(A,I6,A,4I6)') 'setup> for ts ',J1,' +,-,p+,p-:',PLUS(J1),MINUS(J1),POINTERP(J1),POINTERM(J1)
!           ENDDO
!        ENDIF
    
         IF ((.NOT.NOPOINTS).AND.(NATTEMPT.GT.0).AND.(.NOT.DUMMYTST)) THEN
            IF (AMHT) THEN
              WRITE(*,*) 'setup> AVOIDING MOMENT OF INITERIA CALC FOR AMH'
            ELSE
               DO J1=1,NTS
                  READ(UTS,REC=J1) (LOCALPOINTS(J2),J2=1,3*NATOMS)
                  CALL INERTIAWRAPPER(LOCALPOINTS,NATOMS,angleAxis,IXM,IYM,IZM)
!                 IF (DEBUG) WRITE(*,'(A,I6,2F17.7,3I6,3F15.5)') 'setup> ',J1,ETS(J1),FVIBTS(J1),HORDERTS(J1),
!    1                                                            PLUS(J1),MINUS(J1),IXM,IYM,IZM
                  IF ((ABS(IXM-IXTS(J1)).GT.IDIFFTOL).OR.
     1                (ABS(IYM-IYTS(J1)).GT.IDIFFTOL).OR.
     1                (ABS(IZM-IZTS(J1)).GT.IDIFFTOL)) THEN
                     WRITE(*,'(A,I10)') 'setup> WARNING - principal moments of inertia do not agree with input for ts ',J1
                     WRITE(*,'(A)') 'values in ts.data:'
                     WRITE(*,'(3F20.10)') IXTS(J1),IYTS(J1),IZTS(J1)
                     WRITE(*,'(A)') 'values recalculated from points.ts:'
                     WRITE(*,'(3F20.10)') IXM,IYM,IZM
!                    STOP  
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
      ELSE
         WRITE(*,'(A)') 'setup> no transition states found'
         INQUIRE(FILE='ts.data',EXIST=YESNO)
         IF (YESNO) THEN
            PRINT '(A)','ERROR - file ts.data already exists. Will not overwrite.'
            STOP
         ENDIF
         OPEN(UNIT=UTSDATA,FILE='ts.data',STATUS='NEW')
         NTS=0
      ENDIF
!
!  Create a ts entry for DUMMYTS runs if there are minima that seem to lack such entries.
!
      IF (DUMMYTST) THEN
         PRINT '(A)',' setup> shortest distances for local minima:'
         DO J1=1,NMIN
            IF (MINDISTMIN(J1).GT.HUGE(1.0D0)/1.0D1) THEN
               PRINT '(A,I8,A,G20.10)',' setup> in setup, minimum ',J1,' shortest distance unassigned'
               READ(UMIN,REC=J1) (LOCALPOINTS(J2),J2=1,3*NATOMS)
               DO J3=1,NMIN
                  IF (J3.EQ.J1) CYCLE
                  READ(UMIN,REC=J3) (LOCALPOINTS2(J2),J2=1,3*NATOMS)
                  CALL MINPERMDIST(LOCALPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY,
     &                             RMAT,.FALSE.)
                  IF ((DISTANCE.LT.MINDISTMIN(J1)).OR.(DISTANCE.LT.MINDISTMIN(J3))) THEN
                     IF (DISTANCE.LT.MINDISTMIN(J1)) MINDISTMIN(J1)=DISTANCE
                     IF (DISTANCE.LT.MINDISTMIN(J3)) MINDISTMIN(J3)=DISTANCE
!
!  Must create an entry in ts.data in this case.
!  ETS,FVIBTS,HORDERTS,PLUS,MINUS,IXTS,IYTS,IZTS
!
                     IF (IMFRQT) THEN
                        PRINT '()', "setup> ERROR: can''t guess negative eigenvalue - don''t use DUMMYTS and IMFRQ"
                        STOP
                     ENDIF

                     WRITE(UTSDATA,'(2F20.10,3I10,3F20.10)') TSEGUESS(EMIN(J1),EMIN(J3),MINCURVE(J1),MINCURVE(J3),DISTANCE),
     &                           TSFVIBGUESS(EMIN(J1),EMIN(J3),FVIBMIN(J1),FVIBMIN(J3),MINFRQ2(J1),MINFRQ2(J3),NATOMS),
     &                           1,J3,J1,1.0D0,1.0D0,1.0D0
                     CALL FLUSH(UTSDATA,ISTAT)
                     NTS=NTS+1
                     IF (NTS.GT.MAXTS) CALL TSDOUBLE
                     ETS(NTS)=TSEGUESS(EMIN(J1),EMIN(J3),MINCURVE(J1),MINCURVE(J3),DISTANCE)
                     FVIBTS(NTS)=TSFVIBGUESS(EMIN(J1),EMIN(J3),FVIBMIN(J1),FVIBMIN(J3),MINFRQ2(J1),MINFRQ2(J3),NATOMS)
                     HORDERTS(NTS)=1
                     IXTS(NTS)=1.0D0
                     IYTS(NTS)=1.0D0
                     IZTS(NTS)=1.0D0
                     PLUS(NTS)=J3
                     MINUS(NTS)=J1
                     IF (DIJKSTRAT .OR. KSHORTESTPATHST) TSATTEMPT(NTS)=0
                     IF (DEBUG) WRITE(*,'(A,I6,A)') 'setup> dummy ts ',NTS,' writing parameters to file ts.data'
!
!  Update ts pointers.
!
                     POINTERP(NTS)=-1
                     POINTERM(NTS)=-1
                     IF (TOPPOINTER(PLUS(NTS)).GT.0) POINTERP(NTS)=TOPPOINTER(PLUS(NTS))
                     IF (TOPPOINTER(MINUS(NTS)).GT.0) POINTERM(NTS)=TOPPOINTER(MINUS(NTS))
                     TOPPOINTER(PLUS(NTS))=NTS
                     TOPPOINTER(MINUS(NTS))=NTS
                  ENDIF
               ENDDO
            ENDIF
            PRINT '(A,I8,A,G20.10)',' setup> shortest distance for minimum ',J1,' is ',MINDISTMIN(J1)
         ENDDO
      ENDIF
!
!  Set transition state vibrational product to unity for consistency if NOFRQS is true.
!  Won;t work with FREEPAIRS unless the OPTIM runs are also run with NOFRQS.
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
!
!  Procedure to remove selected stationary points specified by min.remove and ts.remove.
!  First line of each file gives the numbers of structures to remove.
!
      IF (REMOVESP) THEN
         OPEN(UNIT=1,FILE='min.remove',STATUS='OLD')
         READ(1,*) NMINREMOVE
         ALLOCATE(MINREMOVE(NMINREMOVE),MINPREV(NMIN))
         READ(1,*) MINREMOVE(1:NMINREMOVE)
         CLOSE(1)
         PRINT '(A)','setup> removing the following minima:'
         PRINT '(10I8)',MINREMOVE(1:NMINREMOVE)
         OPEN(UNIT=2,FILE='min.data.removed',STATUS='UNKNOWN')
         OPEN(UNIT=4,FILE='points.min.removed',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=8*3*NATOMS)
         NDUMMY=0
         MINPREV(1:NMIN)=0
         minloop: DO J1=1,NMIN
            DO J2=1,NMINREMOVE
               IF (MINREMOVE(J2).EQ.J1) THEN
                  PRINT '(A,I8)','setup> removing minimum ',J1
                  CYCLE minloop
               ENDIF
            ENDDO
            PRINT '(A,I8)','setup> not removing minimum ',J1
            NDUMMY=NDUMMY+1
            MINPREV(J1)=NDUMMY
            WRITE(2,'(2F20.10,I6,3F20.10)') EMIN(J1), FVIBMIN(J1), HORDERMIN(J1), IXMIN(J1), IYMIN(J1), IZMIN(J1)
            READ(UMIN,REC=J1) (LOCALPOINTS(J2),J2=1,3*NATOMS)
            WRITE(4,REC=NDUMMY) (LOCALPOINTS(J2),J2=1,3*NATOMS)
         ENDDO minloop
         CLOSE(2)
         CLOSE(4)
!
! rewrite min.A and min.B in min.A.removed and min.B.removed since numbering may change.
!
         NDUMMY=0
         Aloop: DO J1=1,NMINA
            DO J2=1,NMINREMOVE
               IF (MINREMOVE(J2).EQ.LOCATIONA(J1)) THEN
                  PRINT '(A,I8)','setup> removing A minimum ',LOCATIONA(J1)
                  CYCLE Aloop
               ENDIF
            ENDDO
            NDUMMY=NDUMMY+1
         ENDDO Aloop
         OPEN(UNIT=2,FILE='min.A.removed',STATUS='UNKNOWN')
         WRITE(2,'(I8)') NDUMMY
         IF (NDUMMY.EQ.0) THEN
            PRINT '(A)','setup> ERROR - all A minima removed'
            STOP
         ENDIF
         DO J1=1,NMINA
            IF (MINPREV(LOCATIONA(J1)).NE.0) THEN
               WRITE(2,'(I8)') MINPREV(LOCATIONA(J1))
            ENDIF
         ENDDO 
         CLOSE(2)

         NDUMMY=0
         Bloop: DO J1=1,NMINB
            DO J2=1,NMINREMOVE
               IF (MINREMOVE(J2).EQ.LOCATIONB(J1)) THEN
                  PRINT '(A,I8)','setup> removing B minimum ',LOCATIONB(J1)
                  CYCLE Bloop
               ENDIF
            ENDDO
            NDUMMY=NDUMMY+1
         ENDDO Bloop
         OPEN(UNIT=2,FILE='min.B.removed',STATUS='UNKNOWN')
         WRITE(2,'(I8)') NDUMMY
         IF (NDUMMY.EQ.0) THEN
            PRINT '(A)','setup> ERROR - all B minima removed'
            STOP
         ENDIF
         DO J1=1,NMINB
            IF (MINPREV(LOCATIONB(J1)).NE.0) THEN
               WRITE(2,'(I8)') MINPREV(LOCATIONB(J1))
            ENDIF
         ENDDO 
         CLOSE(2)

         OPEN(UNIT=1,FILE='ts.remove',STATUS='OLD')
         READ(1,*) NTSREMOVE
         IF (NTSREMOVE.LE.0) GOTO 444
         ALLOCATE(TSREMOVE(NTSREMOVE))
         READ(1,*) TSREMOVE(1:NTSREMOVE)
         CLOSE(1)
         PRINT '(A)','setup> removing the following transition states:'
         PRINT '(10I8)',TSREMOVE(1:NTSREMOVE)
444      CONTINUE
         OPEN(UNIT=3,FILE='ts.data.removed',STATUS='UNKNOWN')
         OPEN(UNIT=5,FILE='points.ts.removed',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=8*3*NATOMS)
         NDUMMY=0
         tsloop: DO J1=1,NTS
            DO J2=1,NTSREMOVE
               IF (TSREMOVE(J2).EQ.J1) CYCLE tsloop
            ENDDO
            IF (MINPREV(PLUS(J1)).EQ.0) THEN
               PRINT '(A,I8,A,I8,A)','setup> possible ERROR - transition state ',J1,' links minimum ',PLUS(J1), 
     &                               ' which has been removed - removing this ts'
               CYCLE tsloop
            ENDIF
            IF (MINPREV(MINUS(J1)).EQ.0) THEN
               PRINT '(A,I8,A,I8,A)','setup> possible ERROR - transition state ',J1,' links minimum ',MINUS(J1), 
     &                               ' which has been removed - removing this ts'
               CYCLE tsloop
            ENDIF
            NDUMMY=NDUMMY+1
            IF (IMFRQT) THEN
               WRITE(3,'(2F20.10,3I10,4F20.10)') ETS(J1),FVIBTS(J1),HORDERTS(J1),MINPREV(PLUS(J1)),MINPREV(MINUS(J1)),
     &                                        IXTS(J1),IYTS(J1),IZTS(J1),NEGEIG(J1)
            ELSE
               WRITE(3,'(2F20.10,3I10,3F20.10)') ETS(J1),FVIBTS(J1),HORDERTS(J1),MINPREV(PLUS(J1)),MINPREV(MINUS(J1)),
     &                                        IXTS(J1),IYTS(J1),IZTS(J1)
            ENDIF
            READ(UTS,REC=J1) (LOCALPOINTS(J2),J2=1,3*NATOMS)
            WRITE(5,REC=NDUMMY) (LOCALPOINTS(J2),J2=1,3*NATOMS)
         ENDDO tsloop
         CLOSE(3); CLOSE(5)
         STOP
      ENDIF
!
!  Procedure to retain selected stationary points specified by min.retain.
!  First line of each file gives the numbers of structures to retain.
!  All ts linking minima in the retain list are themselves retained.
!
      IF (RETAINSP) THEN
         OPEN(UNIT=1,FILE='min.retain',STATUS='OLD')
         READ(1,*) NMINRETAIN
         ALLOCATE(MINRETAIN(NMINRETAIN),MINPREV(NMIN))
         READ(1,*) MINRETAIN(1:NMINRETAIN)
         CLOSE(1)
         PRINT '(A)','setup> retaining the following minima:'
         PRINT '(10I8)',MINRETAIN(1:NMINRETAIN)
         OPEN(UNIT=2,FILE='min.data.retained',STATUS='UNKNOWN')
         OPEN(UNIT=4,FILE='points.min.retained',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=8*3*NATOMS)
         NDUMMY=0
         MINPREV(1:NMIN)=0
         minloop2: DO J1=1,NMIN
            DO J2=1,NMINRETAIN
               IF (MINRETAIN(J2).EQ.J1) THEN
                  PRINT '(A,I8)','setup> retaining minimum ',J1
                  NDUMMY=NDUMMY+1
                  MINPREV(J1)=NDUMMY
                  WRITE(2,'(2F20.10,I6,3F20.10)') EMIN(J1), FVIBMIN(J1), HORDERMIN(J1), IXMIN(J1), IYMIN(J1), IZMIN(J1)
                  READ(UMIN,REC=J1) (LOCALPOINTS(J3),J3=1,3*NATOMS)
                  WRITE(4,REC=NDUMMY) (LOCALPOINTS(J3),J3=1,3*NATOMS)
                  CYCLE minloop2
               ENDIF
            ENDDO
            PRINT '(A,I8)','setup> removing minimum ',J1
         ENDDO minloop2
         CLOSE(2)
         CLOSE(4)
!
! rewrite min.A and min.B in min.A.retained and min.B.retained since numbering may change.
!
         NDUMMY=0
         Aloop2: DO J1=1,NMINA
            DO J2=1,NMINRETAIN
               IF (MINRETAIN(J2).EQ.LOCATIONA(J1)) THEN
                  NDUMMY=NDUMMY+1
                  PRINT '(A,I8)','setup> retaining A minimum ',LOCATIONA(J1)
                  CYCLE Aloop2
               ENDIF
            ENDDO
         ENDDO Aloop2
         OPEN(UNIT=2,FILE='min.A.retained',STATUS='UNKNOWN')
         WRITE(2,'(I8)') NDUMMY
         IF (NDUMMY.EQ.0) THEN
            PRINT '(A)','setup> ERROR - all A minima removed'
            STOP
         ENDIF
         DO J1=1,NMINA
            IF (MINPREV(LOCATIONA(J1)).NE.0) THEN
               WRITE(2,'(I8)') MINPREV(LOCATIONA(J1))
            ENDIF
         ENDDO 
         CLOSE(2)

         NDUMMY=0
         Bloop2: DO J1=1,NMINB
            DO J2=1,NMINRETAIN
               IF (MINRETAIN(J2).EQ.LOCATIONB(J1)) THEN
                  NDUMMY=NDUMMY+1
                  PRINT '(A,I8)','setup> retaining B minimum ',LOCATIONB(J1)
                  CYCLE Bloop2
               ENDIF
            ENDDO
         ENDDO Bloop2
         OPEN(UNIT=2,FILE='min.B.retained',STATUS='UNKNOWN')
         WRITE(2,'(I8)') NDUMMY
         IF (NDUMMY.EQ.0) THEN
            PRINT '(A)','setup> ERROR - all B minima removed'
            STOP
         ENDIF
         DO J1=1,NMINB
            IF (MINPREV(LOCATIONB(J1)).NE.0) THEN
               WRITE(2,'(I8)') MINPREV(LOCATIONB(J1))
            ENDIF
         ENDDO 
         CLOSE(2)

         OPEN(UNIT=3,FILE='ts.data.retained',STATUS='UNKNOWN')
         OPEN(UNIT=5,FILE='points.ts.retained',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=8*3*NATOMS)
         NDUMMY=0
         tsloop2: DO J1=1,NTS
            DO J2=1,NMINRETAIN
               IF (MINRETAIN(J2).EQ.PLUS(J1)) THEN
                  DO J3=1,NMINRETAIN
                     IF (MINRETAIN(J3).EQ.MINUS(J1)) THEN
                        NDUMMY=NDUMMY+1
                        PRINT '(A,I8,A,2I8)','setup> retaining ts ',J1,' connected to retained minima ',PLUS(J1),MINUS(J1)
                        IF (IMFRQT) THEN
                           WRITE(3,'(2F20.10,3I10,4F20.10)') ETS(J1),FVIBTS(J1),HORDERTS(J1),MINPREV(PLUS(J1)),MINPREV(MINUS(J1)),
     &                                        IXTS(J1),IYTS(J1),IZTS(J1),NEGEIG(J1)
                        ELSE
                           WRITE(3,'(2F20.10,3I10,3F20.10)') ETS(J1),FVIBTS(J1),HORDERTS(J1),MINPREV(PLUS(J1)),MINPREV(MINUS(J1)),
     &                                        IXTS(J1),IYTS(J1),IZTS(J1)
                        ENDIF
                        READ(UTS,REC=J1) (LOCALPOINTS(J4),J4=1,3*NATOMS)
                        WRITE(5,REC=NDUMMY) (LOCALPOINTS(J4),J4=1,3*NATOMS)
                        CYCLE tsloop2
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
            PRINT '(A,I8)','setup> removing ts ',J1
         ENDDO tsloop2
         CLOSE(3); CLOSE(5)
         PRINT '(A,I8,A)','setup> ',NDUMMY,' transition states retained'
         STOP
      ENDIF
!
!  Rate constants.
!
!op226 {{{
!     KMEAN=0.0D0
      IF (ENSEMBLE.EQ.'T') THEN
         DO J1=1,NTS
            KPLUS(J1)  = LOG(1.0D0 * HORDERMIN(PLUS(J1))  / (2.0D0 * PI*HORDERTS(J1))) +
     1             (FVIBMIN(PLUS(J1))  - FVIBTS(J1)) / 2.0D0 - (ETS(J1) - EMIN(PLUS(J1)) )/TEMPERATURE
            IF (FRICTIONT) KPLUS(J1)=KPLUS(J1)+LOG(FRICTIONFAC(NEGEIG(J1)))
            KMINUS(J1) = LOG(1.0D0 * HORDERMIN(MINUS(J1)) / (2.0D0 * PI*HORDERTS(J1))) +
     1             (FVIBMIN(MINUS(J1)) - FVIBTS(J1)) / 2.0D0 - (ETS(J1) - EMIN(MINUS(J1)))/TEMPERATURE
            IF (FRICTIONT) KMINUS(J1)=KMINUS(J1)+LOG(FRICTIONFAC(NEGEIG(J1)))
            IF (ZSYM(1:2).EQ.'CA') KPLUS(J1)=KPLUS(J1)+30.66356D0
            IF (ZSYM(1:2).EQ.'CA') KMINUS(J1)=KMINUS(J1)+30.66356D0
            IF (PLUS(J1).EQ.MINUS(J1)) KPLUS(J1)=KPLUS(J1)+LOG(2.0D0)
            IF (PLUS(J1).EQ.MINUS(J1)) KMINUS(J1)=KMINUS(J1)+LOG(2.0D0)
!           KMEAN=KMEAN+KPLUS(J1)+KMINUS(J1)
            IF (DEBUG) WRITE(*,'(A,3I6,5F15.5,3G20.10)') 'setup> J1,PLUS,MINUS,Ets,E+,E-,k+,k-,<k>=',J1,PLUS(J1),MINUS(J1),
     1                                            ETS(J1),EMIN(PLUS(J1)),EMIN(MINUS(J1)),KPLUS(J1),KMINUS(J1)
         ENDDO
      ELSE
         DO J1=1,NTS
            IF (TOTALE.GT.ETS(J1)) THEN
               KPLUS(J1)  = LOG(1.0D0 * HORDERMIN(PLUS(J1))  / (2*PI*HORDERTS(J1))) +
     1                   (FVIBMIN(PLUS(J1))  - FVIBTS(J1))/2 + (KAPPA-1)*LOG((TOTALE-ETS(J1))/(TOTALE-EMIN(PLUS(J1))))
               KMINUS(J1) = LOG(1.0D0 * HORDERMIN(MINUS(J1)) / (2*PI*HORDERTS(J1))) +
     1                   (FVIBMIN(MINUS(J1)) - FVIBTS(J1))/2 + (KAPPA-1)*LOG((TOTALE-ETS(J1))/(TOTALE-EMIN(MINUS(J1))))
               IF (ZSYM(1:2).EQ.'CA') KPLUS(J1)=KPLUS(J1)+30.66356D0
               IF (ZSYM(1:2).EQ.'CA') KMINUS(J1)=KMINUS(J1)+30.66356D0
               IF (PLUS(J1).EQ.MINUS(J1)) KPLUS(J1)=KPLUS(J1)+LOG(2.0D0)
               IF (PLUS(J1).EQ.MINUS(J1)) KMINUS(J1)=KMINUS(J1)+LOG(2.0D0)
!              KMEAN=KMEAN+KPLUS(J1)+KMINUS(J1)
            ELSE
               KPLUS(J1)=-1.0D250
               KMINUS(J1)=-1.0D250
            ENDIF
         ENDDO
      ENDIF
      !op226 }}}
!     IF (NTS.GT.0) KMEAN=KMEAN/(2.0D0*NTS)
!     PRINT '(A,G20.10)', 'setup> Mean log rate constant=', KMEAN
!
!  Sums of rates out of the intermediate minima
!
!op226 {{{
!       DO J1=1,NMIN
!          KSUM(J1)=0.0D0
!       ENDDO
!       DO J1=1,NTS
!          IF (PLUS(J1).NE.MINUS(J1)) KSUM(PLUS(J1))=KSUM(PLUS(J1))+EXP(KPLUS(J1)-KMEAN)
!          IF (PLUS(J1).NE.MINUS(J1)) KSUM(MINUS(J1))=KSUM(MINUS(J1))+EXP(KMINUS(J1)-KMEAN)
!       ENDDO
!       DO J1=1,NMIN
!          IF (KSUM(J1).GT.0.0D0) THEN
!             KSUM(J1)=LOG(KSUM(J1))+KMEAN
! !           IF (DEBUG) WRITE(*,'(A,I6,2E20.10)') 'setup> J1,KSUM=',J1,KSUM(J1)
!          ENDIF
!       ENDDO
!       DO J1=1,NTS
! !        IF (DEBUG) WRITE(*,'(A,I6,2E20.10)') 'setup> J1,k+,k-=',J1,KPLUS(J1),KMINUS(J1)
!       ENDDO
!op226 }}}
      IF (MERGEDBT) THEN
         CALL MERGEDB
         STOP
      ENDIF
!
!  Add transition states and minima from the <PATHNAME> file.
!  Use GETNEWPATH to do the bookkeeping.
!
!op226 {{{
      IF (ADDPATH) THEN
         CALL MYSYSTEM(STATUS,DEBUG,'cp ' // TRIM(ADJUSTL(PATHNAME)) // ' path.info')
         IF (ADDTRIPLES) THEN
            CALL GETALLPATHS
         ELSE
            CALL GETNEWPATH(0,0)
         ENDIF
      ENDIF
!op226 }}}

      IF (NPFOLD.GT.0) THEN
         INQUIRE(FILE='commit.data',EXIST=YESNO)
         GPFOLD(1:NMIN)=0.0D0
         IF (YESNO) THEN
            PRINT '(A)','setup> Reading initial committor probabilities read from commit.data'
            OPEN(UNIT=1,FILE='commit.data',STATUS='OLD')
            J2=0
            DO J1=1,NMIN
               READ(1,*,END=110) GPFOLD(J1)
               J2=J2+1
            ENDDO 
110         CLOSE(1)
            IF (J2.LT.NMIN) PRINT '(A)','setup> WARNING - end of file commit.data, remaining probabilities set to zero'
         ELSE
            IF (DIRECTION.EQ.'AB') THEN ! GPFOLD is PFA
               DO J1=1,NMINA
                  GPFOLD(LOCATIONA(J1))=1.0D0
               ENDDO
            ELSE ! GPFOLD is PFB
               DO J1=1,NMINB
                  GPFOLD(LOCATIONB(J1))=1.0D0
               ENDDO
            ENDIF
            PRINT '(A)','setup> Initial committor probabilities set to 0 or 1'
!           PRINT '(6G20.10)',GPFOLD(1:NMIN)
         ENDIF
      ENDIF
!
!  Read in the pairs of minima previously searched in pairs.data exists.
!
!{{{
      ALLOCATE(PAIR1(MAXPAIRS),PAIR2(MAXPAIRS))
      INQUIRE(FILE='pairs.data',EXIST=YESNO)
      NPAIRDONE=0
      IF (YESNO) THEN
         OPEN(UNIT=1,FILE='pairs.data',STATUS='OLD')
         DO 
            NPAIRDONE=NPAIRDONE+1
            IF (NPAIRDONE.GT.MAXPAIRS) CALL PAIRDOUBLE
            READ(1,*,END=120) PAIR1(NPAIRDONE), PAIR2(NPAIRDONE)
            IF (DEBUG) PRINT '(A,I8,A,2I8)','setup > previously searched pair number ',
     &                                      NPAIRDONE,' is ',PAIR1(NPAIRDONE), PAIR2(NPAIRDONE)
            IF ((PAIR1(NPAIRDONE).GT.NMIN).OR.(PAIR2(NPAIRDONE).GT.NMIN)) THEN
               PRINT '(A)','setup> ERROR *** minima specified in pairs.data do not exist in min.data'
               STOP
            ENDIF
         ENDDO
120      CLOSE(1)
         NPAIRDONE=NPAIRDONE-1
         PRINT '(A,I8,A)','setup> ',NPAIRDONE,' pairs of minima already searched read from pairs.data'
      ENDIF
! }}}
!
!  Read in the minima previously searched in UNTRAP runs.
!
! {{{
      ALLOCATE(MINDONE(MAXDONE))
      INQUIRE(FILE='min.done',EXIST=YESNO)
      NMINDONE=0
      IF (YESNO) THEN
         OPEN(UNIT=1,FILE='min.done',STATUS='OLD')
         DO 
            NMINDONE=NMINDONE+1
            IF (NMINDONE.GT.MAXDONE) CALL DONEDOUBLE
            READ(1,*,END=121) MINDONE(NMINDONE)
            IF (DEBUG) PRINT '(A,I8,A,2I8)','setup > previously searched minimum ',
     &                                      NMINDONE,' is ',MINDONE(NMINDONE)
         ENDDO
121      CLOSE(1)
         NMINDONE=NMINDONE-1
         PRINT '(A,I8,A)','setup> ',NMINDONE,' minima already searched read from min.done'
      ENDIF
! }}}
!
!  Initialise PAIRDIST array for use in making an intial connection.
!  PAIRDIST should contain zero if the two minima are linked by a transition state.
!
! {{{
      IF (DIJINITT) THEN
         IF (.NOT.INDEXCOSTFUNCTION) THEN
            DO J1=1,NMIN
               READ(UMIN,REC=J1) (LOCALPOINTS(J2),J2=1,3*NATOMS)
               PAIRDIST(J1*(J1+1)/2)=0.0D0
               DO J2=J1+1,NMIN
                  READ(UMIN,REC=J2) (LOCALPOINTS2(J3),J3=1,3*NATOMS)
                  CALL MINPERMDIST(LOCALPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY,
     &                             RMAT,.FALSE.)
                  IF (INTERPCOSTFUNCTION) CALL MINPERMDIST(LOCALPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,
     &                                                     DISTANCE,DIST2,RIGIDBODY,RMAT,INTERPCOSTFUNCTION)
                  PAIRDIST(J2*(J2-1)/2+J1)=DISTANCE
!                 PRINT '(A,3I6,G20.10)','J1,J2,INDEX,DISTANCE=',J1,J2,J2*(J2-1)/2+J1,DISTANCE
               ENDDO
               PRINT '(A,I8)','setup> Finished pair distance calculation for minimum ',J1
               CALL FLUSH(6,ISTAT)
            ENDDO
         ELSE
            DO J1=1,NMIN
               DO J2=J1+1,NMIN
                  PAIRDIST(J2*(J2-1)/2+J1)=1.0D0 ! this is not actually used !
               ENDDO
            ENDDO
         ENDIF
         DO J1=1,NPAIRDONE
            PAIRDIST(MAX(PAIR1(J1),PAIR2(J1))*(MAX(PAIR1(J1),PAIR2(J1))-1)/2+MIN(PAIR1(J1),PAIR2(J1)))=HUGE(1.0D0)
         ENDDO
         DO J1=1,NTS
! JMC n.b. don't apply the nconnmin criteria at this point, hence the huge(1) 's in place of NCONN() for the plus and minus minima.
            CALL CHECKTS(ETS(J1),EMIN(PLUS(J1)),EMIN(MINUS(J1)),KPLUS(J1),KMINUS(J1),HUGE(1),HUGE(1), 
     &                   PLUS(J1),MINUS(J1),.TRUE.,CUT_UNDERFLOW,DEADTS)
            IF (.NOT. DEADTS) THEN
               J2=MAX(PLUS(J1),MINUS(J1))
               J3=MIN(PLUS(J1),MINUS(J1))
               PAIRDIST(J2*(J2-1)/2+J3)=0.0D0
            ENDIF
         ENDDO
!
! Set the pairs for which connections have already been tried to infinite distance,
! so they are not tried again. Don;t overwrite zero distance settings for connections
! that have actually been found!
!
      ENDIF
! }}}
!
! If USEPAIRST is true then read the sequence of minima from file USEPAIRSFILE
! USEPAIRSFILE must be formatted as a single Epath file
!
      IF (USEPAIRST) THEN
         OPEN(UNIT=1,FILE=TRIM(ADJUSTL(USEPAIRSFILE)),STATUS='OLD')
         NUSEPAIRS=0
         DO
            READ(1,*,END=111) NDUMMY, DUMMY, NDUMMY
            NUSEPAIRS=NUSEPAIRS+1
         ENDDO
111      REWIND(1)
         PRINT '(A,A,A,I8,A,I8,A)','setup> Number of lines in file ',TRIM(ADJUSTL(USEPAIRSFILE)),' is ',NUSEPAIRS,' i.e. ',
     &           (NUSEPAIRS+1)/2,' minima'
         NUSEPAIRS=(NUSEPAIRS+1)/2
         ALLOCATE(USEPAIRSMIN(NUSEPAIRS))
         DO J1=1,NUSEPAIRS
            READ(1,*) NDUMMY, DUMMY, USEPAIRSMIN(J1)
            IF (J1.EQ.NUSEPAIRS) EXIT
            READ(1,*) NDUMMY, DUMMY, NDUMMY
         ENDDO
         CLOSE(1)
         PRINT '(A)','setup> Sequence of local minima:'
         PRINT '(15I8)',USEPAIRSMIN(1:NUSEPAIRS)
      ENDIF

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
