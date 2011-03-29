!   CONNECT MODULE IS AN IMPLEMENTATION OF A CONNECTION ALGORITHM FOR FINDING REARRANGEMENT PATHWAYS.
!   COPYRIGHT (C) 2003-2006 SEMEN A. TRYGUBENKO AND DAVID J. WALES
!   THIS FILE IS PART OF CONNECT MODULE. CONNECT MODULE IS PART OF OPTIM.
!
!   OPTIM IS FREE SOFTWARE; YOU CAN REDISTRIBUTE IT AND/OR MODIFY
!   IT UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE AS PUBLISHED BY
!   THE FREE SOFTWARE FOUNDATION; EITHER VERSION 2 OF THE LICENSE, OR
!   (AT YOUR OPTION) ANY LATER VERSION.
!
!   OPTIM IS DISTRIBUTED IN THE HOPE THAT IT WILL BE USEFUL,
!   BUT WITHOUT ANY WARRANTY; WITHOUT EVEN THE IMPLIED WARRANTY OF
!   MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  SEE THE
!   GNU GENERAL PUBLIC LICENSE FOR MORE DETAILS.
!
!   YOU SHOULD HAVE RECEIVED A COPY OF THE GNU GENERAL PUBLIC LICENSE
!   ALONG WITH THIS PROGRAM; IF NOT, WRITE TO THE FREE SOFTWARE
!   FOUNDATION, INC., 59 TEMPLE PLACE, SUITE 330, BOSTON, MA  02111-1307  USA
!
MODULE TRYCONNECTMODULE
     IMPLICIT NONE
     CONTAINS

     SUBROUTINE TRYCONNECT(JS,JF,TSREDO,REDOPATH,REDOPATHXYZ)
          USE NEWNEBMODULE
          USE CONNECTDATA
          USE KEYCONNECT
          USE CONNECTUTILS
          USE NEBTOCONNECT
          USE KEYNEB, ONLY: NIMAGE, NITERMAX
          USE KEY, ONLY: UNRST, DEBUG, FILTH, FILTHSTR, DUMPALLPATHS, TWOD, MAXTSENERGY, RIGIDBODY, &
  &                   PERMDIST, MAXBARRIER, GROWSTRINGT, DESMDEBUG, BULKT
          USE MODGUESS
          USE MODUNRES
          USE MODCHARMM, ONLY: CHRMMT,NCHENCALLS,CHECKOMEGAT,CHECKCHIRALT
          USE MODMEC
          USE KEYUTILS
          USE COMMONS, ONLY: NINTS, PARAM1, PARAM2, PARAM3
          USE PORFUNCS
          USE GSDATA, ONLY: GSITERDENSITY, GSCURITERD=>ITERD
          IMPLICIT NONE
          DOUBLE PRECISION RMAT(3,3)

          INTEGER,INTENT(IN) :: JS,JF
          
          INTEGER         :: I,UNIQUE=0,MINPLUSPOS,MINMINUSPOS,J1
          DOUBLE PRECISION         :: EDUMMY,EDUMMY2,TMPTS(3*NATOMS)
          DOUBLE PRECISION,POINTER :: QPLUS(:),QMINUS(:),EPLUS,EMINUS
          LOGICAL         :: PLUSNEW,MINUSNEW,PATHFAILT,AMIDEFAIL,CHIRALFAIL
          CHARACTER       :: ITSTRING*80, EOFSSTRING*80
          CHARACTER       :: ITSTRINGP*80, EOFSSTRINGP*80
          CHARACTER       :: ITSTRINGM*80, EOFSSTRINGM*80
          DOUBLE PRECISION :: STARTINT(NINTS), FINISHINT(NINTS),DUM
          LOGICAL REDOPATH, REDOPATHXYZ, PERMUTE, EXISTS, NORERUN, PERMONCE
          DOUBLE PRECISION TSREDO(3*NATOMS), QTEMP(3*NATOMS), QLOCAL(3*NATOMS), ETSPREV, ETSDUM, DIST2
          CHARACTER(LEN=5) ZTEMP(NATOMS)
          CHARACTER(LEN=80) DSTRING, TMPSTRING
          INTEGER INVERT, INDEX(NATOMS), J2, IMATCH
          CHARACTER(LEN=2) ZDUM
          CHARACTER(LEN=10) DUMMYS          

          IF (GROWSTRINGT) THEN
             IF (NCONDONE.EQ.1.AND.FCD) THEN
                GSCURITERD = GSITERDENSITY
             ELSE
                GSCURITERD = ITERDENSITY
             ENDIF
          ENDIF

          IF (CHRMMT) NCHENCALLS = 999 ! UPDATE NON-BONDED LIST ON NEXT CALL TO POTENTIAL.

          IF (.NOT.REDOPATH) CALL CHECKPAIR(JS,JF)

          IF (GUESSPATHT) THEN
            IF (UNRST) THEN
               DO J1=1,NRES
                  C(1,J1)=MI(JS)%DATA%X(6*(J1-1)+1)
                  C(2,J1)=MI(JS)%DATA%X(6*(J1-1)+2)
                  C(3,J1)=MI(JS)%DATA%X(6*(J1-1)+3)
                  C(1,J1+NRES)=MI(JS)%DATA%X(6*(J1-1)+4)
                  C(2,J1+NRES)=MI(JS)%DATA%X(6*(J1-1)+5)
                  C(3,J1+NRES)=MI(JS)%DATA%X(6*(J1-1)+6)
               ENDDO
               CALL UPDATEDC
!CALL INT_FROM_CART(.TRUE.,.FALSE.)
!CALL GEOM_TO_VAR(NINTS,STARTINT(1:NINTS)) 
               DO J1=1,NRES
                  C(1,J1)=MI(JF)%DATA%X(6*(J1-1)+1)
                  C(2,J1)=MI(JF)%DATA%X(6*(J1-1)+2)
                  C(3,J1)=MI(JF)%DATA%X(6*(J1-1)+3)
                  C(1,J1+NRES)=MI(JF)%DATA%X(6*(J1-1)+4)
                  C(2,J1+NRES)=MI(JF)%DATA%X(6*(J1-1)+5)
                  C(3,J1+NRES)=MI(JF)%DATA%X(6*(J1-1)+6)
               ENDDO
               CALL UPDATEDC
!CALL INT_FROM_CART(.TRUE.,.FALSE.)
!CALL GEOM_TO_VAR(NINTS,FINISHINT(1:NINTS))
               CALL GUESSPATH(STARTINT,FINISHINT,NINTS,EDIFFTOL,NATOMS)
            ELSE
               CALL GUESSPATH(MI(JS)%DATA%X,MI(JF)%DATA%X,3*NATOMS,EDIFFTOL,NATOMS)
            ENDIF
          ! HOW MANY IMAGES TO USE?
            IF (NIMAGE > IMAGEMAX) PRINT*,'WARNING - NIMAGE IS GREATER THAN IMAGEMAX'
            IF (NIMAGE < 2       ) PRINT*,'WARNING - NIMAGE IS < 2'
            NIMAGE=NINTERP
!           NITERMAX = NIMAGE*ITERDENSITY ! TRY ZERO NEB ITERATIONS IF WE HAVE A GUESSPATH PATH
            NITERMAX = 0
            IF (NINTERP.LT.2) THEN ! NO IMAGES FROM GUESSPATH - REVERT TO USUAL SCHEME
               IF (.NOT.(NCONDONE==1 .AND. FCD)) THEN
                    NIMAGE=IMAGEDENSITY*MI(JF)%DATA%D(JS) &
                        +IMAGEINCR*IMAGEDENSITY*MI(JF)%DATA%D(JS)*(MI(JF)%DATA%NTRIES(JS)-1)
                    IF (NIMAGE >= IMAGEMAX) THEN
                       NIMAGE = IMAGEMAX
!                      MI(JF)%DATA%NTRIES(JS)=NTRIESMAX ! NO POINT TRYING AGAIN WITH THE SAME NUMBER OF IMAGES
                    ENDIF
                    IF (NIMAGE < 2       ) NIMAGE = 2
                    NITERMAX = NIMAGE*ITERDENSITY
               ENDIF

            ENDIF
          ELSEIF (MECCANOT) THEN
          ! HOW MANY IMAGES TO USE?
            NIMAGE=NINT(MIN(MECIMDENS*MI(JF)%DATA%D(JS),MECMAXIMAGES*1.0D0)) ! IMAGE DENSITY TIMES DISTANCE
!           IF (NIMAGE > IMAGEMAX) PRINT*,'WARNING - NIMAGE IS GREATER THAN IMAGEMAX'
            IF (NIMAGE < 1       ) NIMAGE=1
            NITERMAX=NINT(MIN(NIMAGE*MECITDENS,MECMAXIT*1.0D0)) ! NUMBER OF IMAGES TIMES ITERATION DENSITY
            IF (UNRST) THEN
               DO J1=1,NRES
                  C(1,J1)=MI(JS)%DATA%X(6*(J1-1)+1)
                  C(2,J1)=MI(JS)%DATA%X(6*(J1-1)+2)
                  C(3,J1)=MI(JS)%DATA%X(6*(J1-1)+3)
                  C(1,J1+NRES)=MI(JS)%DATA%X(6*(J1-1)+4)
                  C(2,J1+NRES)=MI(JS)%DATA%X(6*(J1-1)+5)
                  C(3,J1+NRES)=MI(JS)%DATA%X(6*(J1-1)+6)
               ENDDO
               CALL UPDATEDC
!CALL INT_FROM_CART(.TRUE.,.FALSE.)
!CALL GEOM_TO_VAR(NINTS,STARTINT(1:NINTS))
               DO J1=1,NRES
                  C(1,J1)=MI(JF)%DATA%X(6*(J1-1)+1)
                  C(2,J1)=MI(JF)%DATA%X(6*(J1-1)+2)
                  C(3,J1)=MI(JF)%DATA%X(6*(J1-1)+3)
                  C(1,J1+NRES)=MI(JF)%DATA%X(6*(J1-1)+4)
                  C(2,J1+NRES)=MI(JF)%DATA%X(6*(J1-1)+5)
                  C(3,J1+NRES)=MI(JF)%DATA%X(6*(J1-1)+6)
               ENDDO
               CALL UPDATEDC
!CALL INT_FROM_CART(.TRUE.,.FALSE.)
!CALL GEOM_TO_VAR(NINTS,FINISHINT(1:NINTS))
               CALL UNMECCANO(.TRUE.,DEBUG,DUM,.FALSE.,STARTINT,FINISHINT,DUM,DUM,DUM,DUM)
            ELSE
               CALL MECCANO(.TRUE.,DEBUG,DUM,.FALSE.,MI(JS)%DATA%X,MI(JF)%DATA%X,DUM,DUM,DUM,DUM)
            ENDIF
            NITERMAX = 0 ! TRY ZERO NEB ITERATIONS IF WE HAVE A MECCANO PATH

          ELSEIF (REDOPATH) THEN
            NIMAGE=1
            NITERMAX = 0 ! TRY ZERO NEB ITERATIONS IF WE HAVE POINTS IN TSREDO
          ELSEIF (NCONDONE==1 .AND. FCD) THEN ! FIRST CYCLE DIFFERENT - PARAMETERS SUPPLIED USING NEWNEB
                                          ! KEYWORD OR NEWNEB DEFAULTS WILL BE USED INSTEAD
               PRINT *, "FIRST CYCLE WILL BE DONE USING EXTERNALLY SUPPLIED PARAMETERS"
          ELSE
            NIMAGE=IMAGEDENSITY*MI(JF)%DATA%D(JS) &
             +IMAGEINCR*IMAGEDENSITY*MI(JF)%DATA%D(JS)*(MI(JF)%DATA%NTRIES(JS)-1)
               IF (NIMAGE >= IMAGEMAX) THEN
                  NIMAGE = IMAGEMAX
!                 MI(JF)%DATA%NTRIES(JS)=NTRIESMAX ! NO POINT TRYING AGAIN WITH THE SAME NUMBER OF IMAGES
               ENDIF
               IF (NIMAGE < 2       ) NIMAGE = 2
               NITERMAX = NIMAGE*ITERDENSITY
          ENDIF
          
          ! BOOK-KEEPING :-)
          IF (.NOT.(MECCANOT.OR.REDOPATH)) THEN             
             IF (GROWSTRINGT) THEN
                WRITE(CHR,'(I7)') INT(NIMAGE*GSCURITERD)
                WRITE(*,'(/1X,A)',ADVANCE='NO') '>>>>>  '//TRIM(ADJUSTL(CHR))//'-ITERATION GS RUN FOR MINIMA '
             ELSE
                WRITE(CHR,'(I7)') NITERMAX
                WRITE(*,'(/1X,A)',ADVANCE='NO') '>>>>>  '//TRIM(ADJUSTL(CHR))//'-ITERATION DNEB RUN FOR MINIMA '
             ENDIF
             WRITE(CHR,'(I5)') JS
             WRITE(*,'(A)',ADVANCE='NO') TRIM(ADJUSTL(CHR))
          
             IF (MI(JS)%DATA%S) THEN
                  WRITE(*,'(A)',ADVANCE='NO') '_S'
             ELSEIF (MI(JS)%DATA%F) THEN
                  WRITE(*,'(A)',ADVANCE='NO') '_F'
             ELSE
                  WRITE(*,'(A)',ADVANCE='NO') '_U'
             ENDIF
             WRITE(CHR,'(I5)') JF
             WRITE(*,'(A)',ADVANCE='NO') ' AND '//TRIM(ADJUSTL(CHR))
             IF (MI(JF)%DATA%S) THEN
                  WRITE(*,'(A)',ADVANCE='NO') '_S'
             ELSEIF (MI(JF)%DATA%F) THEN
                  WRITE(*,'(A)',ADVANCE='NO') '_F'
             ELSE
                  WRITE(*,'(A)',ADVANCE='NO') '_U'
             ENDIF
          
             ! GETTING TS CANDIDATES FROM NEB
             WRITE(CHR,'(I5)') NIMAGE
             WRITE(*,'(A)',ADVANCE='NO') ' USING '//TRIM(ADJUSTL(CHR))//' IMAGES '
             IF (MI(JF)%DATA%NTRIES(JS) > 1) THEN
                  WRITE(CHR,'(I5)') MI(JF)%DATA%NTRIES(JS)
                  WRITE(*,'(A)',ADVANCE='NO') '(ATTEMPT #'//TRIM(ADJUSTL(CHR))//') '
             ENDIF
             WRITE(*,'(A)') ' ...'
          ENDIF

          IF (NIMAGE >= IMAGEMAX) MI(JF)%DATA%NTRIES(JS)=NTRIESMAX ! NO POINT TRYING AGAIN WITH THE SAME NUMBER OF IMAGES
          IF (REDOPATHXYZ) THEN
             NORERUN=.FALSE.
!
! NCONDONE SHOULD BE NTS+1, OTHERWISE ONE PREVIOUS RUN HAS FAILED TO GET A NEW TS FROM THE PATH.<N>.XYZ FILE,
! WHICH SHOULD NEVER HAPPEN?!
!
             IF (NCONDONE.NE.NTS+1) THEN
                PRINT '(2(A,I8))',' TRYCONNECT> ERROR - NCONDONE=',NCONDONE,' NTS=',NTS
                STOP
             ENDIF
             CALL MKFNAMES(NTS+1,FILTH,FILTHSTR,ITSTRING,EOFSSTRING)
             INQUIRE(FILE=TRIM(ADJUSTL(ITSTRING)),EXIST=EXISTS)
             IF (EXISTS) THEN ! ALLOWS FOR RERUN WITH DIFFERENT ENERGY DIFFERENCE CRITERION FOR
                              ! CONSECUTIVE FRAMES IN PATH WITHOUT REDOING ORIGINAL PATH
                PRINT '(2A)',' TRYCONNECT> READING DATA FOR TS FROM EXISTING FILE ',TRIM(ADJUSTL(ITSTRING))
                NORERUN=.TRUE.
             ELSE
                PRINT '(3A)',' TRYCONNECT> NO FILE ',TRIM(ADJUSTL(ITSTRING)),' FOUND'
                REDOPATH=.FALSE.
                REDOPATHXYZ=.FALSE.
                STOP ! REDOPATHXYZ IS ONLY GOING TO WORK AFTER A RUN FROM REDOPOINTS, IN WHICH CASE
                     ! WE SHOULD FINISH AFTER READING THE LAST PATH.
             ENDIF
             IF (NORERUN) THEN
                NTSFOUND=1
                ALLOCATE(TSFOUND(NTSFOUND)%E,TSFOUND(NTSFOUND)%COORD(NOPT),&
   &                     TSFOUND(NTSFOUND)%EVALMIN,TSFOUND(NTSFOUND)%VECS(NOPT))
                OPEN(UNIT=89,FILE=ITSTRING,STATUS='OLD')
                ETSPREV=1.0D100
                TSFOUND(1)%E=-1.0D100
                DO
                   READ(89,*)
                   READ(89,*) DUMMYS, ETSDUM
                   IF ((TSFOUND(1)%E.GT.ETSDUM).AND.(TSFOUND(1)%E.GT.ETSPREV)) GOTO 111
                   DO J1=1,NATOMS
                      READ(89,*) ZDUM,TSFOUND(1)%COORD(3*(J1-1)+1),TSFOUND(1)%COORD(3*(J1-1)+2),TSFOUND(1)%COORD(3*(J1-1)+3)
                   ENDDO
                   ETSPREV=TSFOUND(1)%E
                   TSFOUND(1)%E=ETSDUM
                ENDDO
111             CONTINUE
                IF (DEBUG) PRINT *,'TSFOUND(1)%E=',TSFOUND(1)%E
                PRINT '(A,G20.10)',' TRYCONNECT> ETS=',TSFOUND(1)%E
                CLOSE(89)
             ELSE
                IF (UNRST) THEN 
                   CALL NEWNEB(REDOPATH,TSREDO,-1.0D100,MI(JS)%DATA%X,-1.0D100,MI(JF)%DATA%X,NINTSIN=NINTS)
                ELSE
                   CALL NEWNEB(REDOPATH,TSREDO,-1.0D100,MI(JS)%DATA%X,-1.0D100,MI(JF)%DATA%X)
                ENDIF
             ENDIF
          ELSE
             IF (UNRST) THEN 
                 CALL NEWNEB(REDOPATH,TSREDO,MI(JS)%DATA%E,MI(JS)%DATA%X,MI(JF)%DATA%E,MI(JF)%DATA%X,NINTSIN=NINTS)
             ELSE
                 CALL NEWNEB(REDOPATH,TSREDO,MI(JS)%DATA%E,MI(JS)%DATA%X,MI(JF)%DATA%E,MI(JF)%DATA%X)
             ENDIF
          ENDIF

          ! SAVING NEW TS INTO TS RACK; OTHERWISE - FREE MEMORY IMMEDIATELY
          NTSOLD=NTS
          UNIQUE=0
          DO I=1,NTSFOUND
!              PRINT '(A,2G20.10)',' TRYCONNECT> TSFOUND(I)%E, MAXTSENERGY=',TSFOUND(I)%E, MAXTSENERGY
               AMIDEFAIL=.FALSE.
               IF (CHRMMT .AND. CHECKOMEGAT) &
                  CALL CHECKOMEGA(TSFOUND(I)%COORD,AMIDEFAIL)
               CHIRALFAIL=.FALSE.
               IF (CHRMMT .AND. CHECKCHIRALT) &
                  CALL CHECKCHIRAL(TSFOUND(I)%COORD,CHIRALFAIL)
               IF (CHRMMT .AND. AMIDEFAIL) THEN
                 PRINT '(A,G20.10,A)',' TRYCONNECT> TRANSITION STATE WITH ENERGY ',TSFOUND(I)%E, &
  &                                   ' IGNORED, CIS-TRANS ISOMERISATION OF AN AMIDE BOND DETECTED.'
                 DEALLOCATE(TSFOUND(I)%E,TSFOUND(I)%COORD,TSFOUND(I)%EVALMIN,TSFOUND(I)%VECS)
               ELSEIF (CHRMMT .AND. CHIRALFAIL) THEN
                  PRINT '(A,G20.10,A)',' TRYCONNECT> TRANSITION STATE WITH ENERGY ',TSFOUND(I)%E, &
  &                                    ' IGNORED, INVERSION OF A CHIRAL CA CENTER DETECTED.'
                  DEALLOCATE(TSFOUND(I)%E,TSFOUND(I)%COORD,TSFOUND(I)%EVALMIN,TSFOUND(I)%VECS)
               ELSEIF (TSFOUND(I)%E.GT.MAXTSENERGY) THEN
                  PRINT '(A,G20.10,A)',' TRYCONNECT> TRANSITION STATE WITH ENERGY ',TSFOUND(I)%E,' IGNORED'
                  DEALLOCATE(TSFOUND(I)%E,TSFOUND(I)%COORD,TSFOUND(I)%EVALMIN,TSFOUND(I)%VECS)
               ELSEIF ( ISNEWTS(TSFOUND(I)) ) THEN
                    IF (NTS==TSRACKSIZE) CALL REALLOCATETSRACK
                    NTS=NTS+1; UNIQUE=UNIQUE+1
                    TS(NTS)%DATA%E => TSFOUND(I)%E
                    TS(NTS)%DATA%X => TSFOUND(I)%COORD
                    TS(NTS)%DATA%EVALMIN => TSFOUND(I)%EVALMIN
                    TS(NTS)%DATA%VECS => TSFOUND(I)%VECS
                    TS(NTS)%DATA%BAD=.FALSE.
                    NULLIFY(TSFOUND(I)%E,TSFOUND(I)%COORD,TSFOUND(I)%EVALMIN,TSFOUND(I)%VECS)
               ELSE
                    IF (NCONDONE==1) PRINT *, 'DISCARDED TS #',I
                    DEALLOCATE(TSFOUND(I)%E,TSFOUND(I)%COORD,TSFOUND(I)%EVALMIN,TSFOUND(I)%VECS)
               ENDIF
          ENDDO

          CALL DUMPTS

          ! PRINT INFO AS TO HOW MANY TS ARE ACTUALLY USEFUL
          IF (UNIQUE==NTSFOUND.AND..NOT.UNIQUE==0) THEN
               IF (NTSFOUND==1) THEN
                    WRITE(*,'(1X,A)') 'TS APPEARS TO BE NEW'
               ELSE
                    WRITE(*,'(1X,A)') 'ALL OF TS FOUND APPEAR TO BE NEW'
               ENDIF
          ELSEIF (UNIQUE < NTSFOUND) THEN
               WRITE(CHR,'(I7)') UNIQUE 
               WRITE(*,'(1X,A)') TRIM(ADJUSTL(CHR))//' OF TS FOUND APPEAR TO BE NEW.'
          ELSEIF (UNIQUE ==0 .AND..NOT.NTSFOUND==0) THEN
               WRITE(*,'(1X,A)') 'ALL OF TS FOUND ARE ALREADY KNOWN'
          ENDIF

          ! PATH RUN FOR ALL UNIQUE TS
          DO I=NTS-UNIQUE+1,NTS
               WRITE(CHR,'(I5)') I
               PRINT '(/1X,A)', '>>>>>  PATH RUN FOR TS '//TRIM(ADJUSTL(CHR))//' ...'
               ALLOCATE( QPLUS(NOPT),QMINUS(NOPT),EPLUS,EMINUS )
               CALL MKFNAMES(I,FILTH,FILTHSTR,ITSTRING,EOFSSTRING)
               EDUMMY=TS(I)%DATA%E
               TMPTS=TS(I)%DATA%X
               ! STRUCTURE IN TS(I)%DATA%X IS A STATIONARY POINT WHICH IS WHY WE DON'T NEED TO STORE G AND RMS FOR IT
               GDUMMY(1:3*NATOMS)=0.0D0; RMS=0.0D0 ! WE MUST INITIALIZE THEM HERE, HOWEVER 
               NORERUN=.FALSE.
               IF (REDOPATH) THEN
                  CALL MKFNAMES(NCONDONE,FILTH,FILTHSTR,ITSTRING,EOFSSTRING)
                  IF (REDOPATHXYZ) THEN
                     INQUIRE(FILE=TRIM(ADJUSTL(ITSTRING)),EXIST=EXISTS)
                     IF (EXISTS) THEN ! ALLOWS FOR RERUN WITH DIFFERENT ENERGY DIFFERENCE CRITERION FOR
                                      ! CONSECUTIVE FRAMES IN PATH WITHOUT REDOING ORIGINAL PATH
                        PRINT '(2A)',' TRYCONNECT> READING DATA FOR MINIMA FROM EXISTING FILE ',TRIM(ADJUSTL(ITSTRING))
                        NORERUN=.TRUE.
                     ENDIF
                  ENDIF
               ENDIF
               IF (.NOT.NORERUN) THEN
                  CALL PATH(TMPTS,EDUMMY,GDUMMY,RMS,TS(I)%DATA%EVALMIN,TS(I)%DATA%VECS,  &
                    & .FALSE.,QPLUS,QMINUS,DEBUG,EDUMMY2,EPLUS,EMINUS, &
                    & TS(I)%DATA%SLENGTH,TS(I)%DATA%DISP,TS(I)%DATA%GAMMA,TS(I)%DATA%NTILDE,FRQSTS,FRQSPLUS, &
                    & FRQSMINUS,ITSTRING,EOFSSTRING,PATHFAILT)
               ELSE
                  OPEN(UNIT=89,FILE=ITSTRING,STATUS='OLD')
                  READ(89,*)
                  READ(89,*) DUMMYS, EPLUS
                  DO J1=1,NATOMS
                     READ(89,*) ZDUM,QPLUS(3*(J1-1)+1),QPLUS(3*(J1-1)+2),QPLUS(3*(J1-1)+3)
                  ENDDO
                  DO
                     READ(89,*,END=99)
                     READ(89,*) DUMMYS, EMINUS
                     DO J1=1,NATOMS
                        READ(89,*) ZDUM,QMINUS(3*(J1-1)+1),QMINUS(3*(J1-1)+2),QMINUS(3*(J1-1)+3)
                     ENDDO
                  ENDDO
99                CONTINUE
                  IF (DEBUG) PRINT '(A,G20.10)','EPLUS=',EPLUS
                  IF (DEBUG) PRINT '(A,G20.10)','EMINUS=',EMINUS
                  PRINT '(A,G20.10,A,G20.10)',' TRYCONNECT> E+=',EPLUS,'                      E-=',EMINUS
                  CLOSE(89)
               ENDIF
               DEALLOCATE(TS(I)%DATA%VECS)
               AMIDEFAIL=.FALSE.
               IF (CHRMMT .AND. CHECKOMEGAT) &
                  CALL CHECKOMEGA(QPLUS,AMIDEFAIL)
               IF (CHRMMT .AND. CHECKOMEGAT .AND. .NOT.AMIDEFAIL) &
                  CALL CHECKOMEGA(QMINUS,AMIDEFAIL)
               CHIRALFAIL=.FALSE.
               IF (CHRMMT .AND. CHECKCHIRALT) &
                  CALL CHECKCHIRAL(QPLUS,CHIRALFAIL)
               IF (CHRMMT .AND. CHECKCHIRALT .AND. .NOT.CHIRALFAIL) &
                  CALL CHECKCHIRAL(QMINUS,CHIRALFAIL)
               IF (CHRMMT .AND. AMIDEFAIL) THEN
                  PRINT '(A,G20.10,A)',' TRYCONNECT> TRANSITION STATE WITH ENERGY ',TS(I)%DATA%E, &
  &                                     ' IGNORED, CIS-TRANS ISOMERISATION OF AN AMIDE-BOND DETECTED.'
                  DEALLOCATE(TS(I)%DATA%EVALMIN)
                  DEALLOCATE(QPLUS,QMINUS,EPLUS,EMINUS)
                  TS(I)%DATA%BAD=.TRUE.
                  CYCLE
               ELSEIF (CHRMMT .AND. CHIRALFAIL) THEN
                  PRINT '(A,G20.10,A)',' TRYCONNECT> TRANSITION STATE WITH ENERGY ',TS(I)%DATA%E, &
  &                                        ' IGNORED, INVERSION OF A CHIRAL CA CENTER DETECTED.'
                  DEALLOCATE(TS(I)%DATA%EVALMIN)
                  DEALLOCATE(QPLUS,QMINUS,EPLUS,EMINUS)
                  TS(I)%DATA%BAD=.TRUE.
                  CYCLE
               ELSEIF (PATHFAILT) THEN
                    DEALLOCATE(TS(I)%DATA%EVALMIN)
                    DEALLOCATE(QPLUS,QMINUS,EPLUS,EMINUS)
                    TS(I)%DATA%BAD=.TRUE.
                    CYCLE
               ELSEIF (TS(I)%DATA%E-MAX(EPLUS,EMINUS).GT.MAXBARRIER) THEN
                  PRINT '(2(A,G20.10))',' TRYCONNECT> TRANSITION STATE WITH ENERGY ',TS(I)%DATA%E,' IGNORED, MINIMUM BARRIER=', &
  &                                      TS(I)%DATA%E-MAX(EPLUS,EMINUS)
                  DEALLOCATE(TS(I)%DATA%EVALMIN)
                  DEALLOCATE(QPLUS,QMINUS,EPLUS,EMINUS)
                  TS(I)%DATA%BAD=.TRUE.
                  CYCLE
               ELSE
                  IF (UNRST) CALL TESTSAMEMIN(EPLUS,QPLUS,EMINUS,QMINUS,PATHFAILT)
                  IF (PATHFAILT) THEN
                      DEALLOCATE(TS(I)%DATA%EVALMIN)
                      DEALLOCATE(QPLUS,QMINUS,EPLUS,EMINUS)
                      TS(I)%DATA%BAD=.TRUE.
                      CYCLE 
                  ENDIF
               ENDIF

!              PERMONCE=.FALSE. ! A DEGENERATE REARRANGEMENT OF A PERMUATIONAL ISOMER COULD OTHERWISE
!                               ! CAUSE AN INFINITE LOOP.
333            CALL ISNEWMIN(EPLUS,QPLUS,MINPLUSPOS,PLUSNEW,REDOPATH,PERMUTE,INVERT,INDEX,IMATCH)
!              IF ((.NOT.(PERMUTE.AND.REDOPATH)).OR.PERMONCE.OR.REDOPATHXYZ) THEN
                  CALL ISNEWMIN(EMINUS,QMINUS,MINMINUSPOS,MINUSNEW,REDOPATH,PERMUTE,INVERT,INDEX,IMATCH)
!              ENDIF

!                IF (PERMUTE.AND.REDOPATH.AND.(.NOT.PERMONCE).AND.(.NOT.REDOPATHXYZ)) THEN
!                   PERMONCE=.TRUE.
!                   IF (IMATCH.EQ.2) THEN ! PERMUTE THE FINAL MINIMUM INSTEAD!
!                      PRINT '(A)',' TRYCONNECT> PERMUTING THE FINISH MINIMUM TO ALIGN IT'
!                      DO J2=1,NATOMS
!                         QTEMP(3*(J2-1)+1)=INVERT*MI(IMATCH)%DATA%X(3*(INDEX(J2)-1)+1)
!                         IF (TWOD) THEN
!                            QTEMP(3*(J2-1)+2)=MI(IMATCH)%DATA%X(3*(INDEX(J2)-1)+2)
!                            QTEMP(3*(J2-1)+3)=MI(IMATCH)%DATA%X(3*(INDEX(J2)-1)+3)
!                         ELSE
!                            QTEMP(3*(J2-1)+2)=INVERT*MI(IMATCH)%DATA%X(3*(INDEX(J2)-1)+2)
!                            QTEMP(3*(J2-1)+3)=INVERT*MI(IMATCH)%DATA%X(3*(INDEX(J2)-1)+3)
!                         ENDIF
!                      ENDDO
!                      MI(IMATCH)%DATA%X(1:NOPT)=QTEMP(1:NOPT)
! !
! ! RECALCULATE ALL DISTANCES 
! !
!                      IF (PERMDIST) THEN
!                         CALL MINPERMDIST(MI(1)%DATA%X,MI(2)%DATA%X, NATOMS, &
!   &                                DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,D,DIST2,RIGIDBODY,RMAT)
!                         D=SQRT(D)
!                      ELSE
!                         CALL NEWMINDIST(MI(1)%DATA%X,MI(2)%DATA%X,NATOMS,D,.FALSE.,.FALSE.,'AX   ',.TRUE.,RIGIDBODY,DEBUG,RMAT)
!                      ENDIF
!                      MI(2)%DATA%D(1)=D
!                      IF (INTERPCOSTFUNCTION) MI(2)%DATA%INTERP(1)=INTERPVALUE(MI(1)%DATA%X(:),MI(2)%DATA%X(:),D)
!                      DO J2=3,NMIN
!                         CALL NEWMINDIST(MI(J2)%DATA%X,MI(2)%DATA%X,NATOMS,D,.FALSE.,.FALSE.,'AX   ',.TRUE.,RIGIDBODY,DEBUG,RMAT)
!                         MI(J2)%DATA%D(2)=D
!                         IF (INTERPCOSTFUNCTION) MI(J2)%DATA%INTERP(2)=INTERPVALUE(MI(J2)%DATA%X(:),MI(2)%DATA%X(:),D)
!                      ENDDO
!                      GOTO 333
!                   ELSE
!                      PRINT '(A)',' TRYCONNECT> PERMUTING THE PATHWAY TO ALIGN IT'
!                      DO J2=1,NATOMS
!                         QTEMP(3*(J2-1)+1)=INVERT*TS(I)%DATA%X(3*(INDEX(J2)-1)+1)
!                         IF (TWOD) THEN
!                            QTEMP(3*(J2-1)+2)=TS(I)%DATA%X(3*(INDEX(J2)-1)+2)
!                            QTEMP(3*(J2-1)+3)=TS(I)%DATA%X(3*(INDEX(J2)-1)+3)
!                         ELSE
!                            QTEMP(3*(J2-1)+2)=INVERT*TS(I)%DATA%X(3*(INDEX(J2)-1)+2)
!                            QTEMP(3*(J2-1)+3)=INVERT*TS(I)%DATA%X(3*(INDEX(J2)-1)+3)
!                         ENDIF
!                      ENDDO
!                      TS(I)%DATA%X(1:NOPT)=QTEMP(1:NOPT)
!                      DO J2=1,NATOMS
!                         QTEMP(3*(J2-1)+1)=INVERT*QPLUS(3*(INDEX(J2)-1)+1)
!                         IF (TWOD) THEN
!                            QTEMP(3*(J2-1)+2)=QPLUS(3*(INDEX(J2)-1)+2)
!                            QTEMP(3*(J2-1)+3)=QPLUS(3*(INDEX(J2)-1)+3)
!                         ELSE
!                            QTEMP(3*(J2-1)+2)=INVERT*QPLUS(3*(INDEX(J2)-1)+2)
!                            QTEMP(3*(J2-1)+3)=INVERT*QPLUS(3*(INDEX(J2)-1)+3)
!                         ENDIF
!                      ENDDO
!                      QPLUS(1:NOPT)=QTEMP(1:NOPT)
!                      DO J2=1,NATOMS
!                         QTEMP(3*(J2-1)+1)=INVERT*QMINUS(3*(INDEX(J2)-1)+1)
!                         IF (TWOD) THEN
!                            QTEMP(3*(J2-1)+2)=QMINUS(3*(INDEX(J2)-1)+2)
!                            QTEMP(3*(J2-1)+3)=QMINUS(3*(INDEX(J2)-1)+3)
!                         ELSE
!                            QTEMP(3*(J2-1)+2)=INVERT*QMINUS(3*(INDEX(J2)-1)+2)
!                            QTEMP(3*(J2-1)+3)=INVERT*QMINUS(3*(INDEX(J2)-1)+3)
!                         ENDIF
!                      ENDDO
!                      QMINUS(1:NOPT)=QTEMP(1:NOPT)
! !
! !  WE NEED TO FIX ALL THE POINTS IN THE PATH.N.XYZ FILE AS WELL, OR WE WILL HAVE
! !  PROBLEMS WHEN WE TRY TO RECONSTRUCT THE FULL PATH.
! !
!                      OPEN(UNIT=81,FILE=ITSTRING,STATUS='OLD')
!                      TMPSTRING=TRIM(ADJUSTL("TEMP." // ITSTRING))
!                      OPEN(UNIT=82,FILE=TMPSTRING,STATUS='UNKNOWN')
! 222                  READ(81,'(A80)',END=444) DSTRING
!                      WRITE(82,'(A80)') DSTRING
!                      READ(81,'(A80)') DSTRING
!                      WRITE(82,'(A80)') DSTRING
!                      DO J2=1,NATOMS
!                         READ(81,*) ZTEMP(J2),QTEMP(3*(J2-1)+1),QTEMP(3*(J2-1)+2),QTEMP(3*(J2-1)+3)
!                      ENDDO
!                      DO J2=1,NATOMS
!                         QLOCAL(1)=INVERT*QTEMP(3*(INDEX(J2)-1)+1)
!                         IF (TWOD) THEN
!                            QLOCAL(2)=QTEMP(3*(INDEX(J2)-1)+2)
!                            QLOCAL(3)=QTEMP(3*(INDEX(J2)-1)+3)
!                         ELSE
!                            QLOCAL(2)=INVERT*QTEMP(3*(INDEX(J2)-1)+2)
!                            QLOCAL(3)=INVERT*QTEMP(3*(INDEX(J2)-1)+3)
!                         ENDIF
!                         WRITE(82,'(A2,4X,3F20.10)') ZTEMP(J2),QLOCAL(1),QLOCAL(2),QLOCAL(3)
!                      ENDDO
!                      GOTO 222
! 444                  CLOSE(81)
!                      CLOSE(82)
!                      CALL SYSTEM('MV ' // TMPSTRING // ' ' // ITSTRING)
!                      GOTO 333
!                   ENDIF
!                ENDIF

               EDUMMY=TS(I)%DATA%E
               TMPTS=TS(I)%DATA%X
               IF (DUMPALLPATHS) CALL MAKEALLPATHINFO(TMPTS,QPLUS,QMINUS,EDUMMY,EPLUS,EMINUS,FRQSTS,FRQSPLUS,FRQSMINUS)
               
               WRITE(CHR,'(I7)') MINPLUSPOS
               WRITE(CHR2,'(I7)') MINMINUSPOS
               100 FORMAT (8X,A,T65,A)
               IF ( .NOT.PLUSNEW .AND. .NOT.MINUSNEW ) THEN
                    WRITE(*,100) 'KNOWN (#'//TRIM(ADJUSTL(CHR))//')','KNOWN (#'//TRIM(ADJUSTL(CHR2))//')'
                    CALL NEWCONNECTION(MINPLUSPOS,MINMINUSPOS,I)
                    DEALLOCATE(QPLUS,QMINUS,EPLUS,EMINUS)
                    CALL SETDISTANCE(MINPLUSPOS,MINMINUSPOS,0.0D0)
                    IF (INTERPCOSTFUNCTION) CALL SETINTERP(MINPLUSPOS,MINMINUSPOS,0.0D0)
               ELSE IF ( PLUSNEW .AND. MINUSNEW ) THEN
                    WRITE(CHR2,'(I7)') MINPLUSPOS+1
                    WRITE(*,100) '*NEW* (PLACED IN '//TRIM(ADJUSTL(CHR))//')','*NEW* (PLACED IN '//TRIM(ADJUSTL(CHR2))//')'

                    CALL ADDNEWMIN(EPLUS,QPLUS)
                    CALL ADDNEWMIN(EMINUS,QMINUS)
                    CALL NEWCONNECTION(MINPLUSPOS,MINPLUSPOS+1,I)
                    MI(MINPLUSPOS+1)%DATA%D(MINPLUSPOS)=0.0D0
                    IF (INTERPCOSTFUNCTION) MI(MINPLUSPOS+1)%DATA%INTERP(MINPLUSPOS)=0.0D0
               ELSE IF ( PLUSNEW .OR. MINUSNEW ) THEN
                    IF ( PLUSNEW ) THEN
                         WRITE(*,100) '*NEW* (PLACED IN '//TRIM(ADJUSTL(CHR))//')','KNOWN (#'//TRIM(ADJUSTL(CHR2))//')'
                         CALL ADDNEWMIN(EPLUS,QPLUS)
                         DEALLOCATE(EMINUS,QMINUS)
                    ELSE
                         WRITE(*,100) 'KNOWN (#'//TRIM(ADJUSTL(CHR))//')','*NEW* (PLACED IN '//TRIM(ADJUSTL(CHR2))//')'
                         DEALLOCATE(EPLUS,QPLUS)
                         CALL ADDNEWMIN(EMINUS,QMINUS)
                    ENDIF
                    CALL NEWCONNECTION(MINPLUSPOS,MINMINUSPOS,I)
                    CALL SETDISTANCE(MINPLUSPOS,MINMINUSPOS,0.0D0)
                    IF (INTERPCOSTFUNCTION) CALL SETINTERP(MINPLUSPOS,MINMINUSPOS,0.0D0)
               ENDIF
          ENDDO

     END SUBROUTINE TRYCONNECT

END MODULE TRYCONNECTMODULE
