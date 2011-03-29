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
MODULE NEWCONNECTMODULE
     IMPLICIT NONE
     CONTAINS

!    CALL       NEWCONNECT(NATOMS,EINITIAL,Q,EFINAL,FIN,DIST,.TRUE.,REDOPATH,REDOPATHXYZ)
     SUBROUTINE NEWCONNECT(NA,EII,QQ,EFF,FINFIN,ENDPOINTSEP,PTEST,REDOPATH,REDOPATHXYZ)

          USE CONNECTDATA
          USE KEYCONNECT
          USE CONNECTUTILS
          USE MODCHARMM
          USE DECIDEWHATTOCONNECT
          USE TRYCONNECTMODULE
          USE IDOMODULE
          USE PORFUNCS
          USE AMHGLOBALS, ONLY : NMRES
          USE KEY,ONLY : BHDISTTHRESH, BHINTERPT, BHDEBUG, DIJKSTRALOCAL, DUMPDATAT, REOPTIMISEENDPOINTS, &
  &                      AMHT, SEQ, MIN1REDO, MIN2REDO, REDOE1, REDOE2, BULKT, TWOD, PERMDIST, RIGIDBODY, &
  &                      INTCONSTRAINTT, INTLJT, INTTST
          USE COMMONS,ONLY : PARAM1, PARAM2, PARAM3, ZSYM
          USE MODNEB,ONLY :  NEWCONNECTT
          IMPLICIT NONE

          INTEGER,INTENT(IN)              :: NA
          DOUBLE PRECISION           :: ENDPOINTSEP,EII,EFF,QQ(3*NA),FINFIN(3*NA)
          LOGICAL,INTENT(IN)              :: PTEST
          DOUBLE PRECISION,POINTER              :: EI,EF
          DOUBLE PRECISION,POINTER,DIMENSION(:) :: Q,FIN

          INTEGER :: JS,JF,J1,J2,NSTART,POSITION,J3,NCOUNT
          CHARACTER(LEN=132) :: STR
          LOGICAL REDOPATH, REDOPATHXYZ, YESNO, SUCCESS, MINNEW, PERMUTE, CHANGED, BIGGERGAP, NOPRINT
          DOUBLE PRECISION TSREDO(3*NA), DSTART, DFINISH
          DOUBLE PRECISION, POINTER :: PINTERPCOORDS(:), PENERGY
          DOUBLE PRECISION INTERPCOORDS(3*NA), ENERGY, OLDDISTS, OLDDISTF
          DOUBLE PRECISION CSTART(3*NA), CFINISH(3*NA), ESTART, EFINISH, RMSINITIAL, RMSFINAL, LGDUMMY(3*NA), ETS, RMSTS
          INTEGER, ALLOCATABLE :: TEMPDIJPAIR(:,:)
          DOUBLE PRECISION, ALLOCATABLE :: TEMPDIJPAIRDIST(:)
          INTEGER INVERT, INDEX(NA), IMATCH, NMINSAVE, NMINSAVE2, ISTAT, MYJS, MYJF
          CHARACTER(LEN=80) FNAMEF
          CHARACTER(LEN=20) EFNAME
          DOUBLE PRECISION BHENERGY
          COMMON /BHINTE/ BHENERGY
          DOUBLE PRECISION DIST2, RMAT(3,3)

          FINISHED=.FALSE.
          ALLOCATE(EI,EF,Q(3*NA),FIN(3*NA))
          EI=EII;EF=EFF;Q=QQ;FIN=FINFIN;

          MOREPRINTING=PTEST
          IF (MOREPRINTING) THEN
               CALL ALLKEYCONNECTPRINT
               PRINT*
          ENDIF
          INQUIRE(FILE='REDOPOINTS',EXIST=YESNO)
          IF (YESNO) THEN
             IF (REDOPATH.AND.(.NOT.REDOPATHXYZ)) THEN
                PRINT '(A)',' NEWCONNECT> TRANSITION STATE COORDINATES WILL BE READ FROM FILE REDOPOINTS'
!            ELSE
!               PRINT '(A)',' NEWCONNECT> WARNING - REDOPOINTS FILE PRESENT, BUT NO REDOPATH KEYWORD'
             ENDIF
          ELSE
             IF (REDOPATH.AND.(.NOT.REDOPATHXYZ)) THEN
                PRINT '(A)',' NEWCONNECT> WARNING - REDOPATH KEYWORD WAS SPECIFIED, BUT NO REDOPOINTS FILE'
                REDOPATH=.FALSE.
             ENDIF
          ENDIF
          IF (REDOPATH.AND.(.NOT.REDOPATHXYZ)) THEN
!
!  IT IS IMPORTANT TO USE EII, EFF, QQ AND FINFIN HERE BECAUSE EI, EF, Q AND FIN ARE
!  DECLARED AS POINTERS. TRYING TO ASIGN STATIC DATA TO THE POINTERS CANNOT WORK!!
!
             OPEN(99,FILE='REDOPOINTS',STATUS='OLD')
             ALLOCATE(MIN1REDO(3*NA),MIN2REDO(3*NA))
             IF (AMHT) THEN
                NCOUNT=0
                DO J2=1,NMRES
                   IF (SEQ(J2).EQ.8) THEN
                      NCOUNT=NCOUNT+1
                      READ(99,*) QQ(3*(NCOUNT-1)+1),QQ(3*(NCOUNT-1)+2),QQ(3*(NCOUNT-1)+3)
                      READ(99,*) QQ(3*(NCOUNT-1)+1),QQ(3*(NCOUNT-1)+2),QQ(3*(NCOUNT-1)+3)
                      NCOUNT=NCOUNT+1
                      READ(99,*) QQ(3*(NCOUNT-1)+1),QQ(3*(NCOUNT-1)+2),QQ(3*(NCOUNT-1)+3)
                   ELSE
                      NCOUNT=NCOUNT+1
                      READ(99,*) QQ(3*(NCOUNT-1)+1),QQ(3*(NCOUNT-1)+2),QQ(3*(NCOUNT-1)+3)
                      NCOUNT=NCOUNT+1
                      READ(99,*) QQ(3*(NCOUNT-1)+1),QQ(3*(NCOUNT-1)+2),QQ(3*(NCOUNT-1)+3)
                      NCOUNT=NCOUNT+1
                      READ(99,*) QQ(3*(NCOUNT-1)+1),QQ(3*(NCOUNT-1)+2),QQ(3*(NCOUNT-1)+3)
                   ENDIF
                ENDDO
             ELSE
                READ(99,*) (QQ(J1),J1=1,3*NA)
             ENDIF
             NREDO=1
41           CONTINUE
             IF (AMHT) THEN
                NCOUNT=0
                DO J2=1,NMRES
                   IF (SEQ(J2).EQ.8) THEN
                      NCOUNT=NCOUNT+1
                      READ(99,*,END=42) FINFIN(3*(NCOUNT-1)+1),FINFIN(3*(NCOUNT-1)+2),FINFIN(3*(NCOUNT-1)+3)
                      READ(99,*) FINFIN(3*(NCOUNT-1)+1),FINFIN(3*(NCOUNT-1)+2),FINFIN(3*(NCOUNT-1)+3)
                      NCOUNT=NCOUNT+1
                      READ(99,*) FINFIN(3*(NCOUNT-1)+1),FINFIN(3*(NCOUNT-1)+2),FINFIN(3*(NCOUNT-1)+3)
                   ELSE
                      NCOUNT=NCOUNT+1
                      READ(99,*,END=42) FINFIN(3*(NCOUNT-1)+1),FINFIN(3*(NCOUNT-1)+2),FINFIN(3*(NCOUNT-1)+3)
                      NCOUNT=NCOUNT+1
                      READ(99,*) FINFIN(3*(NCOUNT-1)+1),FINFIN(3*(NCOUNT-1)+2),FINFIN(3*(NCOUNT-1)+3)
                      NCOUNT=NCOUNT+1
                      READ(99,*) FINFIN(3*(NCOUNT-1)+1),FINFIN(3*(NCOUNT-1)+2),FINFIN(3*(NCOUNT-1)+3)
                   ENDIF
                ENDDO
             ELSE
                READ(99,*,END=42) (FINFIN(J1),J1=1,3*NA)
             ENDIF
             NREDO=NREDO+1
             GOTO 41
42           CONTINUE
             CALL POTENTIAL(QQ,EII,LGDUMMY,.TRUE.,.FALSE.,RMSINITIAL,.FALSE.,.FALSE.)
             IF (CHRMMT.AND.ACESOLV) NCHENCALLS=ACEUPSTEP-1
             CALL POTENTIAL(FINFIN,EFF,LGDUMMY,.TRUE.,.FALSE.,RMSFINAL,.FALSE.,.FALSE.)
             EI=EII;EF=EFF;Q=QQ;FIN=FINFIN;
             PRINT '(3(A,I6))',' NEWCONNECT> REDOPOINTS FILE CONTAINS ',(NREDO-1)/2,' TRANSITION STATES AND ',(NREDO+1)/2,' MINIMA'
             PRINT '(A,F20.10)',' NEWCONNECT> START AND FINISH MINIMA OVERWRITTEN FROM REDOPOINTS FILE'
             WRITE(*,'(A,2(G20.10,A))') ' NEWCONNECT> INITIAL ENERGY=',EI,' RMS FORCE=',RMSINITIAL
             WRITE(*,'(A,2(G20.10,A))') ' NEWCONNECT> FINAL   ENERGY=',EF,' RMS FORCE=',RMSFINAL
             CALL MINPERMDIST(QQ,FINFIN,NA,DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,ENDPOINTSEP,DIST2,RIGIDBODY,RMAT)
             WRITE(*,'(A,G20.10)') ' NEWCONNECT> SEPARATION=',ENDPOINTSEP
          ENDIF

          CALL INITIALISE(NA,EI,Q,EF,FIN,ENDPOINTSEP)

          IF (REDOPATHXYZ) PRINT '(A)',' NEWCONNECT> REDO RUN WILL USE AVAILABLE PATH.<N>.XYZ FILES'
          IF (REDOPATH.AND.(.NOT.REDOPATHXYZ)) THEN
             REWIND(99)
             IF (AMHT) THEN
                NCOUNT=0
                DO J2=1,NMRES
                   IF (SEQ(J2).EQ.8) THEN
                      NCOUNT=NCOUNT+1
                      READ(99,*) MIN2REDO(3*(NCOUNT-1)+1),MIN2REDO(3*(NCOUNT-1)+2),MIN2REDO(3*(NCOUNT-1)+3)
                      READ(99,*) MIN2REDO(3*(NCOUNT-1)+1),MIN2REDO(3*(NCOUNT-1)+2),MIN2REDO(3*(NCOUNT-1)+3)
                      NCOUNT=NCOUNT+1
                      READ(99,*) MIN2REDO(3*(NCOUNT-1)+1),MIN2REDO(3*(NCOUNT-1)+2),MIN2REDO(3*(NCOUNT-1)+3)
                   ELSE
                      NCOUNT=NCOUNT+1
                      READ(99,*) MIN2REDO(3*(NCOUNT-1)+1),MIN2REDO(3*(NCOUNT-1)+2),MIN2REDO(3*(NCOUNT-1)+3)
                      NCOUNT=NCOUNT+1
                      READ(99,*) MIN2REDO(3*(NCOUNT-1)+1),MIN2REDO(3*(NCOUNT-1)+2),MIN2REDO(3*(NCOUNT-1)+3)
                      NCOUNT=NCOUNT+1
                      READ(99,*) MIN2REDO(3*(NCOUNT-1)+1),MIN2REDO(3*(NCOUNT-1)+2),MIN2REDO(3*(NCOUNT-1)+3)
                   ENDIF
                ENDDO
             ELSE
                READ(99,*) (MIN2REDO(J1),J1=1,3*NA)
             ENDIF
          ENDIF

          DO NCONDONE=1,NCONMAX
               WRITE(CHR,'(I5)') NCONDONE
               WRITE(STR,'(3A)') '>>>>>>>>>>>>>>>>>>>>> CONNECT CYCLE ',TRIM(ADJUSTL(CHR)),' >>>>>>>>>>>>>>>>>>>>>'
               WRITE(CHR,'(I5)') NMIN
               WRITE(STR,'(A)')  TRIM(STR)//' '//TRIM(ADJUSTL(CHR))//' MINIMA AND'
               WRITE(CHR,'(I5)') NTS
               WRITE(STR,'(A)')  TRIM(STR)//' '//TRIM(ADJUSTL(CHR))//' TS ARE KNOWN'
               WRITE(*,'(1X,A)') TRIM(STR)//' '//REPEAT('>',107-LEN_TRIM(STR))
!
!  READ IN TRANSITION STATE COORDINATES FROM REDOPOINTS IF AVAILABLE.
!  MUST NOT USE THE POINTERS EI AND EF HERE!!!
!
               IF (REDOPATH.AND.(.NOT.REDOPATHXYZ)) THEN
                  IF (AMHT) THEN
                     MIN1REDO(1:3*NA)=MIN2REDO(1:3*NA)
                     NCOUNT=0
                     DO J2=1,NMRES
                        IF (SEQ(J2).EQ.8) THEN
                           NCOUNT=NCOUNT+1
                           READ(99,*,END=32) TSREDO(3*(NCOUNT-1)+1),TSREDO(3*(NCOUNT-1)+2),TSREDO(3*(NCOUNT-1)+3)
                           READ(99,*) TSREDO(3*(NCOUNT-1)+1),TSREDO(3*(NCOUNT-1)+2),TSREDO(3*(NCOUNT-1)+3)
                           NCOUNT=NCOUNT+1
                           READ(99,*) TSREDO(3*(NCOUNT-1)+1),TSREDO(3*(NCOUNT-1)+2),TSREDO(3*(NCOUNT-1)+3)
                        ELSE
                           NCOUNT=NCOUNT+1
                           READ(99,*,END=32) TSREDO(3*(NCOUNT-1)+1),TSREDO(3*(NCOUNT-1)+2),TSREDO(3*(NCOUNT-1)+3)
                           NCOUNT=NCOUNT+1
                           READ(99,*) TSREDO(3*(NCOUNT-1)+1),TSREDO(3*(NCOUNT-1)+2),TSREDO(3*(NCOUNT-1)+3)
                           NCOUNT=NCOUNT+1
                           READ(99,*) TSREDO(3*(NCOUNT-1)+1),TSREDO(3*(NCOUNT-1)+2),TSREDO(3*(NCOUNT-1)+3)
                        ENDIF
                     ENDDO
                     NCOUNT=0
                     DO J2=1,NMRES
                        IF (SEQ(J2).EQ.8) THEN
                           NCOUNT=NCOUNT+1
                           READ(99,*) MIN2REDO(3*(NCOUNT-1)+1),MIN2REDO(3*(NCOUNT-1)+2),MIN2REDO(3*(NCOUNT-1)+3)
                           READ(99,*) MIN2REDO(3*(NCOUNT-1)+1),MIN2REDO(3*(NCOUNT-1)+2),MIN2REDO(3*(NCOUNT-1)+3)
                           NCOUNT=NCOUNT+1
                           READ(99,*) MIN2REDO(3*(NCOUNT-1)+1),MIN2REDO(3*(NCOUNT-1)+2),MIN2REDO(3*(NCOUNT-1)+3)
                        ELSE
                           NCOUNT=NCOUNT+1
                           READ(99,*) MIN2REDO(3*(NCOUNT-1)+1),MIN2REDO(3*(NCOUNT-1)+2),MIN2REDO(3*(NCOUNT-1)+3)
                           NCOUNT=NCOUNT+1
                           READ(99,*) MIN2REDO(3*(NCOUNT-1)+1),MIN2REDO(3*(NCOUNT-1)+2),MIN2REDO(3*(NCOUNT-1)+3)
                           NCOUNT=NCOUNT+1
                           READ(99,*) MIN2REDO(3*(NCOUNT-1)+1),MIN2REDO(3*(NCOUNT-1)+2),MIN2REDO(3*(NCOUNT-1)+3)
                        ENDIF
                     ENDDO
                     IF (DEBUG) PRINT '(A,I6,A)',' NEWCONNECT> ',NCOUNT,' TS COORDINATES READ FROM REDOPOINTS'
                     CALL POTENTIAL(MIN1REDO,EII,LGDUMMY,.TRUE.,.FALSE.,RMSINITIAL,.FALSE.,.FALSE.)
                     CALL POTENTIAL(TSREDO,ETS,LGDUMMY,.TRUE.,.FALSE.,RMSTS,.FALSE.,.FALSE.)
                     CALL POTENTIAL(MIN2REDO,EFF,LGDUMMY,.TRUE.,.FALSE.,RMSFINAL,.FALSE.,.FALSE.)
                     WRITE(*,'(A,2(G20.10,A))') ' NEWCONNECT> INITIAL ENERGY=',EII,' RMS FORCE=',RMSINITIAL
                     WRITE(*,'(A,2(G20.10,A))') ' NEWCONNECT> TS ENERGY     =',ETS,' RMS FORCE=',RMSTS
                     WRITE(*,'(A,2(G20.10,A))') ' NEWCONNECT> FINAL ENERGY  =',EFF,' RMS FORCE=',RMSFINAL
                     REDOE1=EII
                     REDOE2=EFF
                  ELSE
                     MIN1REDO(1:3*NA)=MIN2REDO(1:3*NA)
                     READ(99,*,END=32) (TSREDO(J1),J1=1,3*NA)
                     READ(99,*) (MIN2REDO(J1),J1=1,3*NA)
                     CALL POTENTIAL(MIN1REDO,EII,LGDUMMY,.TRUE.,.FALSE.,RMSINITIAL,.FALSE.,.FALSE.)
                     CALL POTENTIAL(TSREDO,ETS,LGDUMMY,.TRUE.,.FALSE.,RMSTS,.FALSE.,.FALSE.)
                     IF (CHRMMT.AND.ACESOLV) NCHENCALLS=ACEUPSTEP-1
                     CALL POTENTIAL(MIN2REDO,EFF,LGDUMMY,.TRUE.,.FALSE.,RMSFINAL,.FALSE.,.FALSE.)
                     WRITE(*,'(A,2(G20.10,A))') ' NEWCONNECT> INITIAL ENERGY=',EII,' RMS FORCE=',RMSINITIAL
                     WRITE(*,'(A,2(G20.10,A))') ' NEWCONNECT> TS ENERGY     =',ETS,' RMS FORCE=',RMSTS
                     WRITE(*,'(A,2(G20.10,A))') ' NEWCONNECT> FINAL ENERGY  =',EFF,' RMS FORCE=',RMSFINAL
                     REDOE1=EII
                     REDOE2=EFF
                  ENDIF
                  WRITE(*,'(A)') ' NEWCONNECT> MINIMUM AND TRANSITION STATE COORDINATES READ FROM FILE REDOPOINTS'
                  GOTO 33
32                WRITE(*,'(A)') ' NEWCONNECT> NO MORE COORDINATES IN REDOPOINTS'
                  REDOPATH=.FALSE.
                  CLOSE(99)
33                CONTINUE
               ENDIF

               IF (REDOPATH) THEN
                  CALL TRYCONNECT(1,2,TSREDO,REDOPATH,REDOPATHXYZ,.FALSE.,.FALSE.)
               ELSE
                  CALL DECIDE

!    CALL TRYCONNECT FOR EACH PAIR OF MINIMA SPECIFIED IN THE SHORTEST PATH RETURNED BY DIJKSTRA
!    CALL BHINTERP TO FILL IN ADDTIONAL MINIMA IF REQUIRED.

                  IF (BHINTERPT) THEN
                     NSTART=1
                     CHANGED=.FALSE.
34                   CONTINUE
                     NOPRINT=.FALSE.
                     DO J1=NSTART,NDIJPAIRS
                        IF(.NOT.NOPRINT) PRINT '(A)',' NEWCONNECT> DIJPAIR DISTANCES AND MINIMA:'
                        IF(.NOT.NOPRINT) PRINT '(I8,G20.10,I8,G20.10,I8,G20.10)',(J2,DIJPAIRDIST(J2),DIJPAIR(J2,1), &
                                         MI(DIJPAIR(J2,1))%DATA%E,DIJPAIR(J2,2),MI(DIJPAIR(J2,2))%DATA%E,J2=1,NDIJPAIRS)
                        DO J2=1,NDIJPAIRS
                           IF (DIJPAIRDIST(J2).LT.GEOMDIFFTOL) THEN
!                             PRINT '(A,I8,2G20.10)','J2,DIJPAIRDIST(J2),GEOMDIFFTOL=',J2,DIJPAIRDIST(J2),GEOMDIFFTOL
                              PRINT '(A,G20.10,A,G20.10)',' NEWCONNECT> WARNING - DISTANCE ',DIJPAIRDIST(J2), &
  &                                       ' IS LESS THAN GEOMDIFFTOL=',GEOMDIFFTOL 
!                             STOP
                           ENDIF
                        ENDDO
!
!  ONLY TRY INTERPOLATION AT THE FIRST CONNECTION ATTEMPT FOR EACH PAIR.
!
                        IF ((BHDEBUG).AND.(.NOT.NOPRINT)) PRINT '(A,I8,A,I8,A,I8)', &
  &                      ' BHINTERP> TRIES FOR DIJPAIR MIN ',DIJPAIR(J1,1),' AND MIN ',DIJPAIR(J1,2),' = ', &
  &                        MI(MAX(DIJPAIR(J1,1),DIJPAIR(J1,2)))%DATA%NTRIES(MIN(DIJPAIR(J1,1),DIJPAIR(J1,2)))
                        IF ((MI(MAX(DIJPAIR(J1,1),DIJPAIR(J1,2)))%DATA%NTRIES(MIN(DIJPAIR(J1,1),DIJPAIR(J1,2))).EQ.0).AND. &
                            (DIJPAIRDIST(J1).GT.BHDISTTHRESH)) THEN
                           PRINT '(A,I8,A,I8,A,I8,A,G20.10)',' NEWCONNECT> GAP ',J1,' TRYING BH INTERPOLATION FOR MINIMA   ', &
  &                                DIJPAIR(J1,1),' AND ',DIJPAIR(J1,2),' DIST=',DIJPAIRDIST(J1)
!
!  COUNT THIS AS A CONNECTION ATTEMPT. OTHERWISE WE MAY TRY TO INTERPOLATE TWICE.
!
                           MI(MAX(DIJPAIR(J1,1),DIJPAIR(J1,2)))%DATA%NTRIES(MIN(DIJPAIR(J1,1),DIJPAIR(J1,2)))= &
  &                          MI(MAX(DIJPAIR(J1,1),DIJPAIR(J1,2)))%DATA%NTRIES(MIN(DIJPAIR(J1,1),DIJPAIR(J1,2)))+1
                           CSTART(1:3*NA)=MI(DIJPAIR(J1,1))%DATA%X(1:3*NA)
                           CFINISH(1:3*NA)=MI(DIJPAIR(J1,2))%DATA%X(1:3*NA)
                           ESTART=MI(DIJPAIR(J1,1))%DATA%E
                           EFINISH=MI(DIJPAIR(J1,2))%DATA%E
                           CALL BHINTERP(CSTART,CFINISH,3*NA,NA,INTERPCOORDS,SUCCESS,DSTART,DFINISH,ENERGY,ESTART,EFINISH, &
  &                                      DIJPAIRDIST(J1))
                           BIGGERGAP=.FALSE.
! BS360                    IF ((DSTART.LT.(1.5D0*DIJPAIRDIST(J1))).AND.(DFINISH.LT.(1.5D0*DIJPAIRDIST(J1)))) THEN
                           IF ((DSTART.LT.DIJPAIRDIST(J1)).OR.(DFINISH.LT.DIJPAIRDIST(J1))) THEN
                              IF (SUCCESS.AND.BHDEBUG) PRINT '(A,2G20.10,A)', &
  &                              ' NEWCONNECT> AT LEAST BH ONE INTERPOLATED DISTANCE IS SHORTER: ',DSTART,DFINISH,' ACCEPT'
!
! IF WE ARE GOING TO ALLOW THIS, THEN WE SKIP BISECTION OF EITHER NEW GAP.
!
                               IF ((DSTART.GT.DIJPAIRDIST(J1)).OR.(DFINISH.GT.DIJPAIRDIST(J1))) BIGGERGAP=.TRUE.
                           ELSE IF ((DSTART.GT.DIJPAIRDIST(J1)).AND.(DFINISH.GT.DIJPAIRDIST(J1))) THEN
                              IF (SUCCESS.AND.BHDEBUG) PRINT '(A,2G20.10,A)', &
  &                              ' NEWCONNECT> BH BOTH INTERPOLATED DISTANCES ARE LARGER: ',DSTART,DFINISH,' REJECT'
                              SUCCESS=.FALSE.
                           ENDIF
                           IF (SUCCESS) THEN
!
! IF THE EXTRA MINIMUM IS NEW ADD IT TO THE STACK ALONG WITH ALL ITS DISTANCES.
! MUST ALSO ADD IT TO THE DIJPAIR LIST.
!
                              NULLIFY(PINTERPCOORDS,PENERGY)
                              ALLOCATE(PINTERPCOORDS(3*NA),PENERGY)
!
! HERE WE ARE READING DATA INTO THE POINTERS PENERGY AND PINTERPCOORDS SO THAT WHEN
! THEY ARE SUBSEQUENTLY NULLIFIED THE DATA IN MI%DATA IS PRESERVED. 
! THIS IS A TRULY HORRIBLE HACK CAUSED BY THE INAPPROPRIATE POINTER ATTRIBUTE FOR STATIC DATA.
!
                              OPEN(UNIT=781,FILE='BHSCRATCH',STATUS='UNKNOWN')
                              WRITE(781,*) ENERGY,INTERPCOORDS
                              REWIND(781)
                              READ(781,*) PENERGY,PINTERPCOORDS
                              CLOSE(781)
                              CALL ISNEWMIN(PENERGY,PINTERPCOORDS,POSITION,MINNEW,REDOPATH,PERMUTE,INVERT,INDEX,IMATCH)
                              IF (MINNEW) THEN
                                 WRITE(*,'(A,I7)') ' NEWCONNECT> INTERPOLATED MINIMUM IS NEW MINIMUM ',POSITION
                                 CALL ADDNEWMIN(PENERGY,PINTERPCOORDS)
!                                DO J2=1,NMIN-1
!                                   PRINT '(A,2I8,G20.10)','MIN1 MIN2 D=',NMIN,J2,MI(NMIN)%DATA%D(J2)
!                                ENDDO
                              ELSE
                                 WRITE(*,'(A,I7)') ' NEWCONNECT> INTERPOLATED MINIMUM IS OLD MINIMUM ',POSITION
                                  
                                 IF (DIJPAIR(J1,1).GE.POSITION) THEN
                                    IF (MI(DIJPAIR(J1,1))%DATA%NTRIES(POSITION).GE.NTRIESMAX) THEN
                                       PRINT '(A,I8,A,I8,A,I8,A,G20.10)', &
  &                                      ' BHINTERP> MAXIMUM CONNECTION ATTEMPTS FOR ', &
  &                                      DIJPAIR(J1,1),' AND ',POSITION,' IS ', &
  &                                      MI(DIJPAIR(J1,1))%DATA%NTRIES(POSITION), &
  &                                      ' SKIP; DISTANCE=',MI(DIJPAIR(J1,1))%DATA%D(POSITION)
                                       GOTO 2
                                    ENDIF
                                 ELSE
                                    IF (MI(POSITION)%DATA%NTRIES(DIJPAIR(J1,1)).GE.NTRIESMAX) THEN
                                       PRINT '(A,I8,A,I8,A,I8,A,G20.10)', &
  &                                      ' NEWCONNECT> MAXIMUM CONNECTION ATTEMPTS FOR ', &
  &                                      DIJPAIR(J1,1),' AND ',POSITION,' IS ', &
  &                                      MI(POSITION)%DATA%NTRIES(DIJPAIR(J1,1)), &
  &                                      ' SKIP; DISTANCE=',MI(POSITION)%DATA%D(DIJPAIR(J1,1))
                                       GOTO 2
                                    ENDIF
                                 ENDIF
                                 IF (DIJPAIR(J1,2).GE.POSITION) THEN
                                    IF (MI(DIJPAIR(J1,2))%DATA%NTRIES(POSITION).GE.NTRIESMAX) THEN
                                       PRINT '(A,I8,A,I8,A,I8,A,G20.10)', &
  &                                      ' NEWCONNECT> MAXIMUM CONNECTION ATTEMPTS FOR ', &
  &                                      DIJPAIR(J1,2),' AND ',POSITION,' IS ', &
  &                                      MI(DIJPAIR(J1,2))%DATA%NTRIES(POSITION), &
  &                                      ' SKIP; DISTANCE=',MI(DIJPAIR(J1,2))%DATA%D(POSITION)
                                       GOTO 2
                                    ENDIF
                                 ELSE
                                    IF (MI(POSITION)%DATA%NTRIES(DIJPAIR(J1,2)).GE.NTRIESMAX) THEN
                                       PRINT '(A,I8,A,I8,A,I8,A,G20.10)', &
  &                                      ' NEWCONNECT> MAXIMUM CONNECTION ATTEMPTS FOR ', &
  &                                      DIJPAIR(J1,2),' AND ',POSITION,' IS ', &
  &                                      MI(POSITION)%DATA%NTRIES(DIJPAIR(J1,2)), &
  &                                      ' SKIP; DISTANCE=',MI(POSITION)%DATA%D(DIJPAIR(J1,2))
                                       GOTO 2
                                    ENDIF
                                 ENDIF
!
!  REJECT INTERPOLATED MINIMUM IF IT IS ALREADY IN THE LIST, TO AVOID AN
!  INFINITE LOOP.
!
                                 DO J3=1,NDIJPAIRS
                                    IF ((POSITION.EQ.DIJPAIR(J3,1)).OR.(POSITION.EQ.DIJPAIR(J3,2))) THEN
                                       PRINT '(A)',' NEWCONNECT> INTERPOLATED MINIMUM IS ALREADY IN THE LIST - REJECT'
                                       GOTO 2
                                    ENDIF
                                 ENDDO
                              ENDIF

                              CHANGED=.TRUE.

                              NULLIFY(PINTERPCOORDS,PENERGY) 
!
!  IF THE INTERPOLATED MINMUM IS OLD THERE ARE FOUR POSIBILITIES. IT COULD ALREADY
!  BE CONNECTED TO BOTH ENDPOINTS, TO ONE OF THEM, OR TO NEITHER. THE ENDPOINTS SHOULD
!  NOT HAVE BEEN CHOSEN IN THE FIRST PLACE IF A CONNECTION EXISTS VIA THE INTERPOLATED OLD
!  MINIMUM, SO THAT LEAVES THREE POSSIBILITIES.
!
                              IF (.NOT.MINNEW) THEN
                                 IF (POSITION.GT.DIJPAIR(J1,1)) THEN
                                    OLDDISTS=MI(POSITION)%DATA%D(DIJPAIR(J1,1))
                                 ELSE
                                    OLDDISTS=MI(DIJPAIR(J1,1))%DATA%D(POSITION)
                                 ENDIF
                                 IF (POSITION.GT.DIJPAIR(J1,2)) THEN
                                    OLDDISTF=MI(POSITION)%DATA%D(DIJPAIR(J1,2))
                                 ELSE
                                    OLDDISTF=MI(DIJPAIR(J1,2))%DATA%D(POSITION)
                                 ENDIF
                                 PRINT '(2(A,I8),2(A,G20.10))','NEWCONNECT> OLD DISTANCE BETWEEN ',POSITION,' AND ',DIJPAIR(J1,1), &
  &                                                       ' IS ',OLDDISTS,' NEW DIST=',DSTART
                                 PRINT '(2(A,I8),2(A,G20.10))','NEWCONNECT> OLD DISTANCE BETWEEN ',POSITION,' AND ',DIJPAIR(J1,2), &
  &                                                       ' IS ',OLDDISTF,' NEW DIST=',DFINISH
                                 DSTART=OLDDISTS
                                 DFINISH=OLDDISTF
                              ENDIF

                              IF (DSTART.LT.GEOMDIFFTOL) THEN ! REPLACE ENTRY BY OLD MINIMUM/SECOND MINIMUM PAIR
                                 DIJPAIR(J1,1)=POSITION
                                 DIJPAIRDIST(J1)=DFINISH
                                 NSTART=J1
                                 IF (BIGGERGAP) NSTART=J1+1
                                 GOTO 34 ! GO BACK
                              ELSE IF (DFINISH.LT.GEOMDIFFTOL) THEN ! REPLACE ENTRY BY FIRST MINIMUM/OLD MINIMUM PAIR
                                 DIJPAIR(J1,2)=POSITION
                                 DIJPAIRDIST(J1)=DSTART
                                 NSTART=J1
                                 IF (BIGGERGAP) NSTART=J1+1
                                 GOTO 34 ! GO BACK
                              ELSE ! ADD AN ENTRY TO THE LIST - NEITHER MINIMUM CONNECTED TO THE MINIMUM IN QUESTION
                                 ALLOCATE(TEMPDIJPAIR(NDIJPAIRS,2),TEMPDIJPAIRDIST(NDIJPAIRS))
                                 TEMPDIJPAIR(1:NDIJPAIRS,1:2)=DIJPAIR(1:NDIJPAIRS,1:2)
                                 TEMPDIJPAIRDIST(1:NDIJPAIRS)=DIJPAIRDIST(1:NDIJPAIRS)
                                 DEALLOCATE(DIJPAIR,DIJPAIRDIST)
                                 NDIJPAIRS=NDIJPAIRS+1
                                 ALLOCATE(DIJPAIR(NDIJPAIRS,2),DIJPAIRDIST(NDIJPAIRS))
                                 DO J2=1,J1-1
                                    DIJPAIR(J2,1:2)=TEMPDIJPAIR(J2,1:2)
                                    DIJPAIRDIST(J2)=TEMPDIJPAIRDIST(J2)
                                 ENDDO
                                 DIJPAIR(J1,1)=TEMPDIJPAIR(J1,1)
                                 DIJPAIR(J1,2)=POSITION
                                 DIJPAIR(J1+1,1)=POSITION
                                 DIJPAIR(J1+1,2)=TEMPDIJPAIR(J1,2)
                                 DIJPAIRDIST(J1)=DSTART
                                 DIJPAIRDIST(J1+1)=DFINISH
                                 DO J2=J1+1,NDIJPAIRS-1
                                    DIJPAIR(J2+1,1:2)=TEMPDIJPAIR(J2,1:2)
                                    DIJPAIRDIST(J2+1)=TEMPDIJPAIRDIST(J2)
                                 ENDDO
                                 NSTART=J1
                                 IF (BIGGERGAP) NSTART=J1+2
                                 DEALLOCATE(TEMPDIJPAIR,TEMPDIJPAIRDIST)
                                 GOTO 34 ! GO BACK, SINCE NDIJPAIRS HAS INCREASED BY ONE
                              ENDIF
                           ELSE
                              NSTART=J1+1 ! TO AVOID AN INFINITE LOOP
                           ENDIF
                        ELSE
                           NOPRINT=.TRUE.
                        ENDIF
2                       CONTINUE
                     ENDDO
                     IF (CHANGED) THEN
                        WRITE(*,'(A)') ' NEWCONNECT> THE UNCONNECTED MINIMA IN THE CHAIN AND THEIR DISTANCES ARE NOW:'
                        DO J1=1,NDIJPAIRS
                           WRITE(*,'(I6,F12.2,I6)',ADVANCE="NO") DIJPAIR(J1,1),DIJPAIRDIST(J1),DIJPAIR(J1,2)
                        ENDDO
                     ENDIF
                  ENDIF
                  WRITE(*,'(A)') ' ' ! TO ADVANCE TO THE NEXT LINE

                  IF (.NOT.NEWCONNECTT) THEN ! THIS MUST BE A BHINTERP ONLY RUN
                     PRINT '(A,I8)',' NEWCONNECT> NMIN=',NMIN
                     IF (DUMPDATAT) THEN ! CALL GEOPT TO GENERATE MIN.DATA.INFO ENTRIES
                        FNAMEF='POINTS.FINAL.BH'
                        EFNAME='  '
                        DO J1=3,NMIN ! DON'T INCLUDE THE END POINT MINIMA
!                          WRITE(*,*) 'IN NEWCONNECT.F90, BHENERGY=',BHENERGY
                           CALL FLUSH(6,ISTAT)
                           REOPTIMISEENDPOINTS=.FALSE.
                           BHENERGY=MI(J1)%DATA%E
                           CALL GEOPT(FNAMEF,EFNAME,MI(J1)%DATA%X(1:3*NA)) ! FILE NAMES ARE DUMMIES HERE
                           CALL FLUSH(881,ISTAT)    ! FLUSH MIN.DATA.INFO UNIT AGAIN
                        ENDDO
                     ENDIF
                     STOP
                  ENDIF
!
! IF INTCONSTRAINT IS SPECIFIED USE THE INTLBFGS ROUTINE TO GENERATE A SERIES
! OF IMAGES, MINIMISE THEM, AND USE THE MINIMA IN ORDER FOR SEPARATE DNEB RUNS WITH THE
! TRUE POTENTIAL.
! OR, USE THE INTERPOLATION ROUTINE TO GENERATE TRANSITION STATES DIRECTLY.
!
                  DO J1=1,NDIJPAIRS
                     NMINSAVE=NMIN
                     JS=MAX(DIJPAIR(J1,1),DIJPAIR(J1,2))
                     JF=MIN(DIJPAIR(J1,1),DIJPAIR(J1,2))
!
! CHECK THAT WE HAVEN'T ALREADY CONNECTED THIS PAIR IN A PREVIOUS SEARCH IN THIS CYCLE.
!
                     IF (DIJPAIRDIST(J1).EQ.0.0D0) THEN
                        PRINT '(A,I8,A,2I8,A,G20.10)',' NEWCONNECT> SKIPPING DIJKSTRA PAIR NUMBER ',J1,' MINIMA ', &
  &                            JS,JF,' DISTANCE IN NOW',DIJPAIRDIST(J1)
                        CYCLE
                     ENDIF
                     IF (MI(JS)%DATA%D(JF).EQ.0.0D0) THEN
                        PRINT '(A,I8,A,2I8,A,G20.10)',' NEWCONNECT> SKIPPING DIJKSTRA PAIR NUMBER ',J1,' MINIMA ', &
  &                            JS,JF,' DISTANCE IS NOW ',MI(JS)%DATA%D(JF)
                        CYCLE
                     ENDIF
!                    WRITE(*,'(A,I5,A,2I5)') ' NEWCONNECT> TRYING DIJKSTRA PAIR NUMBER ',J1,' MINIMA ',JS,JF
!
! DON'T COUNT CONSTRAINED POTENTIAL INTERPOLATIONS AS CONNECTION ATTEMPTS.
!
                     IF (.NOT.(INTCONSTRAINTT.OR.INTLJT)) THEN
                        MI(JS)%DATA%NTRIES(JF)=MI(JS)%DATA%NTRIES(JF)+1
                     ELSEIF (.NOT.(INTERPCOSTFUNCTION)) THEN
                        MI(JS)%DATA%NTRIES(JF)=MI(JS)%DATA%NTRIES(JF)+1
                     ELSE
!
!    SET INTERP WEIGHT TO INFINITY TO PREVENT RETRIES.
!
                        MI(JS)%DATA%INTERP(JF)=HUGE(1.0D0)
                     ENDIF
!
! MINIMUM JF BECOMES IMAGE NUMBER ONE AND MINIMUM JS BECOMES IMAGE NIMAGE+2
! JUST TO BE CONFUSING, TRYCONNECT HAS THE ARGUMENTS AS JS,JF, SO THE MINIMA
! GET SWAPPED AROUND IN MAKING THE DNEB INTERPOLATION. THIS SHOULDN'T MATTER,
! EXCEPT THAT THE RECONNECTION PROCEDURE FOR INTCONSTRAINTT/INTLJT BELOW NEEDS TO
! KNOW THIS ORDER!
!
                     CALL TRYCONNECT(JF,JS,TSREDO,REDOPATH,REDOPATHXYZ,INTCONSTRAINTT,INTLJT)
                     IF ((MI(JS)%DATA%NTRIES(JF).GE.NTRIESMAX).AND.(MI(JS)%DATA%D(JF).GT.1.0D-5)) THEN
                        MI(JS)%DATA%D(JF)=HUGE(1.0D0)
                        IF (INTERPCOSTFUNCTION) MI(JS)%DATA%INTERP(JF)=HUGE(1.0D0)
                     ENDIF

!                    PRINT '(A,I6,A,I6)','NEWCONNECT> AFTER TRYCONNECT FOR MINIMA ',JS,' AND ',JF
!                    PRINT '(A,G20.10)','NEWCONNECT> D=',MI(JS)%DATA%D(JF)
!                    PRINT '(A,G20.10)','NEWCONNECT> INTERP=',MI(JS)%DATA%INTERP(JF)
!
! REWEIGHT DISTANCE METRIC BETWEEN MINIMA JS, JF AND ALL THE NEW MINIMA FOUND IN THE
! MOST RECENT CONNECTION ATTEMPT. THIS IS TO ENCOURAGE TRYING NEW CONNECTION ATTEMPTS
! WITHIN THIS SET OF MINIMA, EVEN IF THE DISTANCES WOULD OTHERWISE CAUSE SOME OTHER PAIRS
! TO BE TRIED.
!
                     ! PRINT '(A,2I8)',' NEWCONNECT> NMINSAVE,NMIN=',NMINSAVE,NMIN
                     IF (DIJKSTRALOCAL.NE.1.0D0) THEN
                        DO J2=NMINSAVE+1,NMIN
                           MI(J2)%DATA%D(JS)=MI(J2)%DATA%D(JS)*DIJKSTRALOCAL
                           MI(J2)%DATA%D(JF)=MI(J2)%DATA%D(JF)*DIJKSTRALOCAL
                           PRINT '(A,3I8,2G20.10)','NEWCONNECT> J2,JS,JF=',J2,JS,JF
                           PRINT '(A,3I8,2G20.10)','NEWCONNECT> J2,JS,JF,D(JS),D(JF)=',J2,JS,JF,MI(J2)%DATA%D(JS),MI(J2)%DATA%D(JF)
                           DO J3=J2+1,NMIN
                              MI(J3)%DATA%D(J2)=MI(J3)%DATA%D(J2)*DIJKSTRALOCAL
                              PRINT '(A,2I8,G20.10)','NEWCONNECT> J3,J2=',J2,J3,MI(J3)%DATA%D(J2)
                           ENDDO
                        ENDDO
                     ENDIF
!                    PRINT '(A)','NEWCONNECT DISTANCES SUMMARY:'
!                    DO J2=1,NMIN
!                       DO J3=J2+1,NMIN
!                          PRINT '(A,2I8,G20.10)','NEWCONNECT> J3,J2=',J2,J3,MI(J3)%DATA%D(J2)
!                       ENDDO
!                    ENDDO
!
! DO PROPER DNEB SEARCHES FOR NEW MINIMA FOUND IN THE PREVIOUS RUN IN ORDER.
!
                     IF ((INTCONSTRAINTT.OR.INTLJT).AND.(.NOT.INTTST)) THEN
                        FCD=.FALSE. ! TURN OFF SPECIAL FIRST CONNECTION IMAGE SPECIFICATION IF SET
                        NMINSAVE2=NMIN ! BECAUSE NMIN OCULD INCREASE IN THE NEXT TRYCONNECT CALLS! 
                        IF (NMINSAVE2-NMINSAVE.GT.0) THEN
                           MYJS=MAX(JF,NMINSAVE+1)
                           MYJF=MIN(JF,NMINSAVE+1)
!
! DON'T TRY THE PAIR IF THEY HAVE ALREADY BEEN CONNECTED BY A PREVIOUS RUN!
!
                           MI(MYJS)%DATA%NTRIES(MYJF) = MI(MYJS)%DATA%NTRIES(MYJF) + 1
                           IF (MI(MYJS)%DATA%D(MYJF).GT.0.0D0) THEN
                              WRITE(*,'(A,2I6)') ' NEWCONNECT> TRYING MINIMA FROM INTERPOLATION RUN: ',MYJS,MYJF
                              CALL TRYCONNECT(MYJF,MYJS,TSREDO,REDOPATH,REDOPATHXYZ,.FALSE.,.FALSE.)
                              IF ((MI(MYJS)%DATA%NTRIES(MYJF).GE.NTRIESMAX).AND.  &
                                  (MI(MYJS)%DATA%D(MYJF).GT.1.0D-5)) MI(MYJS)%DATA%D(MYJF)=HUGE(MI(MYJS)%DATA%D(MYJF))
                           ENDIF
                           DO J2=1,NMINSAVE2-NMINSAVE-1
                              MYJF=NMINSAVE+J2
                              MYJS=NMINSAVE+J2+1
                              IF (MI(MYJS)%DATA%D(MYJF).GT.0.0D0) THEN
                                 WRITE(*,'(A,2I6)') ' NEWCONNECT> TRYING MINIMA FROM INTERPOLATION RUN: ',MYJS,MYJF
                                 MI(MYJS)%DATA%NTRIES(MYJF) = MI(MYJS)%DATA%NTRIES(MYJF) + 1
!
!    SET EDGE WEIGHT TO INFINITY IF WE HAVE REACHED THE MAXIMUM NUMBER OF TRIES FOR THIS PAIR.
!
                                CALL TRYCONNECT(MYJF,MYJS,TSREDO,REDOPATH,REDOPATHXYZ,.FALSE.,.FALSE.)
                                IF ((MI(MYJS)%DATA%NTRIES(MYJF).GE.NTRIESMAX).AND.  &
                                    (MI(MYJS)%DATA%D(MYJF).GT.1.0D-5)) MI(MYJS)%DATA%D(MYJF) = HUGE(MI(MYJS)%DATA%D(MYJF))
                              ENDIF
                           ENDDO
                           MYJS=NMINSAVE2
                           MYJF=JS
                           IF (MI(MYJS)%DATA%D(MYJF).GT.0.0D0) THEN
                              WRITE(*,'(A,2I6)') ' NEWCONNECT> TRYING MINIMA FROM INTERPOLATION RUN: ',MYJS,MYJF
                              MI(MYJS)%DATA%NTRIES(MYJF) = MI(MYJS)%DATA%NTRIES(MYJF) + 1
                              CALL TRYCONNECT(MYJF,MYJS,TSREDO,REDOPATH,REDOPATHXYZ,.FALSE.,.FALSE.)
                              IF ((MI(MYJS)%DATA%NTRIES(MYJF).GE.NTRIESMAX).AND.  &
                                  (MI(MYJS)%DATA%D(MYJF).GT.1.0D-5)) MI(MYJS)%DATA%D(MYJF) = HUGE(MI(MYJS)%DATA%D(MYJF))
                           ENDIF
                        ELSE
                           MYJS=JS
                           MYJF=JF
                           IF (MI(MYJS)%DATA%D(MYJF).GT.0.0D0) THEN
                              WRITE(*,'(A,2I6)') ' NEWCONNECT> TRYING MINIMA FROM CONSTRAINED POTENTIAL RUN: ',MYJS,MYJF
                              MI(MYJS)%DATA%NTRIES(MYJF) = MI(MYJS)%DATA%NTRIES(MYJF) + 1
                              CALL TRYCONNECT(MYJF,MYJS,TSREDO,REDOPATH,REDOPATHXYZ,.FALSE.,.FALSE.)
                              IF ((MI(MYJS)%DATA%NTRIES(MYJF).GE.NTRIESMAX).AND.  &
                                  (MI(MYJS)%DATA%D(MYJF).GT.1.0D-5)) MI(MYJS)%DATA%D(MYJF) = HUGE(MI(MYJS)%DATA%D(MYJF))
                           ENDIF
                        ENDIF
                        FCD=.TRUE.
                     ENDIF
                  ENDDO
                  DEALLOCATE(DIJPAIR,DIJPAIRDIST)
               ENDIF
               IF (FINISHED) EXIT
          ENDDO

          CALL OUTPUT
          CALL DEINITIALISE
          DEALLOCATE(Q,FIN,EI,EF)
          IF (ALLOCATED(MIN1REDO)) DEALLOCATE(MIN1REDO)
          IF (ALLOCATED(MIN2REDO)) DEALLOCATE(MIN2REDO)
     END SUBROUTINE NEWCONNECT

END MODULE NEWCONNECTMODULE
