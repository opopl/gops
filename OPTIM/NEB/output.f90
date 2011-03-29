!   NEB MODULE IS AN IMPLEMENTATION OF THE NUDGED ELASTIC BAND METHOD FOR PERFORMING DOUBLE-ENDED PATHWAY SEARCHES.
!   COPYRIGHT (C) 2003-2006 SEMEN A. TRYGUBENKO AND DAVID J. WALES
!   THIS FILE IS PART OF NEB MODULE. NEB MODULE IS PART OF OPTIM.
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
MODULE NEBOUTPUT
     IMPLICIT NONE
     CONTAINS

     SUBROUTINE TSLOCATOR
          USE KEY,ONLY: BFGSTST,UNRST,NSTEPS,MACHINE, GROWSTRINGT, INTEPSILON, REDOTSIM
          USE GSDATA, ONLY: EVOLVESTRINGT
          USE KEYOUTPUT
          USE MODCHARMM
          USE NEBDATA
          USE KEYNEB,ONLY:NIMAGE,DEBUG
          USE NEBTOCONNECT
          USE CHARUTILS
          USE MODGUESS
          USE MODMEC
          USE LINKEDLIST
          USE MODEFOL
          USE INTCOMMONS, ONLY : DESMINT, NINTC, NNZ, KD, INTNEWT
          USE COMMONS, ONLY : REDOPATH, REDOPATHNEB

          IMPLICIT NONE
          
          INTEGER :: I,J,NT,ITDONE=0,J1,RECLEN
          INTEGER,PARAMETER :: MAXPRINTOUT = 50, ITMAX  = 30
          DOUBLE PRECISION :: EDUMMY,EVALMIN,EVALMAX,MAXE,VECSNORM
          LOGICAL :: TSCONVERGED,T
          DOUBLE PRECISION,DIMENSION(3*NATOMS) :: LGDUMMY, VECS, DIAG
          INTEGER :: MLOC
          DOUBLE PRECISION :: TIME, TIME0
          DOUBLE PRECISION :: DPRAND
          LOGICAL :: KNOWE, KNOWG, KNOWH ! JMC
          COMMON /KNOWN/ KNOWE, KNOWG, KNOWH ! JMC
          CHARACTER(LEN=256) :: FILENAME, METHSTR
          TYPE(CHAIN),POINTER :: FIRST,DUMMY

          LOGICAL :: TMPINTNEWT, FAILED

          NULLIFY(FIRST,DUMMY)
          NT = 0

          IF (REDOPATHNEB) THEN
             NT=1
             MAXE=-1.0D100
             MLOC=REDOTSIM+1
             PRINT '(A,F20.10)',' TSLOCATOR> TRANSITION STATE HAS ENERGY ',EEE(REDOTSIM+1)
             ALLOCATE(FIRST)
             DUMMY=>FIRST 
             DUMMY%I=MLOC
             NULLIFY(DUMMY%NEXT)
          ELSE
             SELECT CASE(CANDIDATES) ! IDENTIFY TS CANDIDATES
             CASE('HIGH')
                  NT=1
                  MAXE=-1.0D100
                  PRINT *,'NIMAGE=',NIMAGE
                  DO J1=2,NIMAGE+1
                     PRINT '(A,I6,2G20.10)','OUTPUT> J1,EEE,MAXE=',J1,EEE(J1),MAXE
                     IF (EEE(J1).GT.MAXE) THEN
                        MLOC=J1
                        MAXE=EEE(J1)
                     ENDIF
                  ENDDO
                  ALLOCATE(FIRST)
                  DUMMY=>FIRST 
                  DUMMY%I=MLOC
                  NULLIFY(DUMMY%NEXT)
             CASE('ALL','MAXIM')
                  DO I=2,NIMAGE+1 
                       T=.TRUE.
                       IF (CANDIDATES=='MAXIM') THEN
                            IF ( EEE(I-1) < EEE(I) .AND. EEE(I) > EEE(I+1) ) THEN
                                 T=.TRUE.
                            ELSE
                                 T=.FALSE.
                            ENDIF
                       ENDIF
!                      PRINT '(A,I6,3F20.10,L5)','I,EEE(I-1),EEE(I),EEE(I+1),T=',I,EEE(I-1),EEE(I),EEE(I+1),T
                       IF (T) THEN
                            NT=NT+1
                            IF (ASSOCIATED(FIRST)) THEN
                                 ALLOCATE(DUMMY%NEXT)
                                 DUMMY=>DUMMY%NEXT
                            ELSE
                                 ALLOCATE(FIRST)
                                 DUMMY=>FIRST
                            ENDIF
                            DUMMY%I = I ! IS A POSITION OF A MAXIMUM IN ARRAY XYZ
                            NULLIFY(DUMMY%NEXT)
                       ENDIF
                  ENDDO
             END SELECT
          ENDIF

          IF (ASSOCIATED(FIRST)) THEN
               DUMMY=>FIRST
          ELSE
               ! THIS SHOULD COPE WITH HIGHLY ASYMMETRIC PROFILES, WHICH OTHERWISE LOOKS MONOTONIC
               ! UNTIL WE TRY A HUGE NUMBER OF IMAGES. 
               PRINT '(1X,A)', 'NO MAXIMUM IN PROFILE - USING HIGHEST IMAGE'
               ALLOCATE(FIRST)
               DUMMY=>FIRST
               IF (EEE(2).GT.EEE(NIMAGE+1)) THEN
                  DUMMY%I = 2
               ELSE
                  DUMMY%I = NIMAGE+1
               ENDIF
               NULLIFY(DUMMY%NEXT)
          ENDIF

!         WRITE(*,'(1X,A)',ADVANCE='NO') 'FOLLOWING IMAGES ARE CANDIDATES FOR TS:'
!         DO J=1,NTSMAX
!              IF (J<MAXPRINTOUT) WRITE(*,'(I5)',ADVANCE='NO') DUMMY%I-1
!              IF (.NOT.ASSOCIATED(DUMMY%NEXT)) EXIT
!              DUMMY=>DUMMY%NEXT
!         ENDDO
!         WRITE(*,'(A)') '.'

          ! ------ BS360 : MORE GENERAL PRINTOUT (WITHOUT MAXPRINTOUT) -------
          WRITE(*,'(1X,A,I4,A)',ADVANCE='NO') 'FOLLOWING ',NT,' IMAGES ARE CANDIDATES FOR TS:'
          DO J=1,NTSMAX
               WRITE(*,'(I5)',ADVANCE='NO') DUMMY%I-1
               IF (.NOT.ASSOCIATED(DUMMY%NEXT)) EXIT
               DUMMY=>DUMMY%NEXT
               !MSB50
          ENDDO
          PRINT *,' '
          ! ------ END BS360 ---------------------------
          
          IF (OPTIMIZETS) THEN
               WRITE(*,'(1X,A)',ADVANCE='NO') 'CONVERGED TO TS (NUMBER OF ITERATIONS):     '

               DUMMY=>FIRST
               NTSFOUND=0
               CALL MYCPU_TIME(STARTTIME,.FALSE.)
               DO J=1,NTSMAX
                    CALL MYCPU_TIME(TIME0,.FALSE.)
                    EDUMMY=EEE(DUMMY%I)
                    LGDUMMY(1:3*NATOMS)=TRUEGRAD((DUMMY%I-1)*3*NATOMS+1:DUMMY%I*3*NATOMS)
                    KNOWE=.TRUE.
                    KNOWG=.TRUE.
                    IF (REDOPATH) THEN
                       KNOWG = .FALSE.
                       KNOWE = .FALSE.
                    ENDIF
                    IF (BFGSTST) THEN
                       IF (UNRST) THEN ! JMC
                          KNOWG=.FALSE. ! IS THIS NEEDED NOW THAT GDUMMY IS SET? DJW
                          VECS(1:NINTS)=TANVEC(1:NINTS,DUMMY%I-1)
                          VECSNORM=SUM(VECS(1:NINTS)**2)
                          IF (VECSNORM.EQ.0.0D0) THEN  ! JUST IN CASE TANVEC IS SOMEHOW NOT SET? E.G. FOR REDOPATH !
                             IF (DEBUG) PRINT '(A)', ' OUTPUT> SETTING RANDOM INITIAL VECTOR FOR EIGENVECTOR'
                             DO J1=1,NINTS
                                VECS(J1)=DPRAND()*2-1.0D0
                             ENDDO
                             CALL VECNORM(VECS,NINTS)
                          ENDIF
                          CALL INTBFGSTS(NSTEPS,XYZ(NOPT*(DUMMY%I-1)+1:NOPT*DUMMY%I),  &
                   &       EDUMMY,LGDUMMY,TSCONVERGED,RMS,EVALMIN,EVALMAX,VECS,ITDONE,.TRUE.,DEBUG)
                       ELSE
                          IF (DESMINT) THEN
                             TMPINTNEWT = INTNEWT
                             INTNEWT = .FALSE. ! LINEAR TRANSFORMATION ONLY
                             ! CONVERT INTERNAL TANGENTS TO CARTESIANS
                             CALL TRANSBACKDELTA(TANVEC(1:NOPT,DUMMY%I-1),VECS,XYZCART(3*NATOMS*(DUMMY%I-1)+1:3*NATOMS*DUMMY%I), &
                                  & NINTC,3*NATOMS,NNZ,KD,FAILED,DEBUG,INTEPSILON)                             
                             INTNEWT = TMPINTNEWT
                             VECSNORM=SUM(VECS(1:3*NATOMS)**2)
                             IF (VECSNORM.EQ.0.0D0) THEN  ! TANVEC ISN'T SET FOR GUESSPATH, MECCANO, UNMECCANO
                                IF (DEBUG) PRINT '(A)', ' OUTPUT> SETTING RANDOM INITIAL VECTOR FOR EIGENVECTOR'
                                DO J1=1,3*NATOMS
                                   VECS(J1)=DPRAND()*2-1.0D0
                                ENDDO
                                CALL VECNORM(VECS,3*NATOMS)
                             ENDIF
                          ELSE
                             VECS(1:NOPT)=TANVEC(1:NOPT,DUMMY%I-1)
                             VECSNORM=SUM(VECS(1:NOPT)**2)
                             IF (VECSNORM.EQ.0.0D0) THEN  ! TANVEC ISN'T SET FOR GUESSPATH, MECCANO, UNMECCANO
                                IF (DEBUG) PRINT '(A)', ' OUTPUT> SETTING RANDOM INITIAL VECTOR FOR EIGENVECTOR'
                                DO J1=1,NOPT
                                   VECS(J1)=DPRAND()*2-1.0D0
                                ENDDO
                                CALL VECNORM(VECS,NOPT)
                             ENDIF
                          ENDIF
                          IF (GROWSTRINGT.OR.REDOPATH) THEN
                             KNOWG = .FALSE.
                             KNOWE = .FALSE.
                          ENDIF

                          IF (DESMINT) THEN
                             CALL BFGSTS(NSTEPS,XYZCART(3*NATOMS*(DUMMY%I-1)+1:3*NATOMS*DUMMY%I),  &
                                  &       EDUMMY,LGDUMMY,TSCONVERGED,RMS,EVALMIN,EVALMAX,VECS,ITDONE,.TRUE.,PRINTOPTIMIZETS)
                          ELSE
                             CALL BFGSTS(NSTEPS,XYZ(NOPT*(DUMMY%I-1)+1:NOPT*DUMMY%I),  &
                                  &       EDUMMY,LGDUMMY,TSCONVERGED,RMS,EVALMIN,EVALMAX,VECS,ITDONE,.TRUE.,PRINTOPTIMIZETS)
                          ENDIF
                       ENDIF
                    ELSE
                       IF (DESMINT) THEN
                          CALL EFOL(XYZCART(3*NATOMS*(DUMMY%I-1)+1:3*NATOMS*DUMMY%I),TSCONVERGED, &
                               &   NSTEPS,EDUMMY,ITDONE,EVALMIN,DEBUG,DIAG,2)
                       ELSE
                          CALL EFOL(XYZ(NOPT*(DUMMY%I-1)+1:NOPT*DUMMY%I),TSCONVERGED, &
                               &   NSTEPS,EDUMMY,ITDONE,EVALMIN,DEBUG,DIAG,2)
                       ENDIF
                    ENDIF
                    CALL MYCPU_TIME(TIME,.FALSE.)

!                   IF (CHRMMT) CALL CHECKTS(DUMMY,EVALMIN,TSCONVERGED) ! THIS IS NOW A DUMMY ROUTINE!

                    IF (TSCONVERGED) THEN
                         NTSFOUND=NTSFOUND+1
                         IF (DESMINT) THEN
                            ALLOCATE(TSFOUND(NTSFOUND)%E,TSFOUND(NTSFOUND)%COORD(3*NATOMS),&
                                 &TSFOUND(NTSFOUND)%EVALMIN,TSFOUND(NTSFOUND)%VECS(3*NATOMS))
                            TSFOUND(NTSFOUND)%VECS=VECS(1:3*NATOMS)
                            TSFOUND(NTSFOUND)%COORD=XYZCART(3*NATOMS*(DUMMY%I-1)+1:3*NATOMS*DUMMY%I)
                         ELSE
                            ALLOCATE(TSFOUND(NTSFOUND)%E,TSFOUND(NTSFOUND)%COORD(NOPT),&
                                 &TSFOUND(NTSFOUND)%EVALMIN,TSFOUND(NTSFOUND)%VECS(NOPT))
                            TSFOUND(NTSFOUND)%VECS=VECS(1:NOPT)
                            TSFOUND(NTSFOUND)%COORD=XYZ(NOPT*(DUMMY%I-1)+1:NOPT*DUMMY%I)
                         ENDIF
                         TSFOUND(NTSFOUND)%E=EDUMMY
                         TSFOUND(NTSFOUND)%EVALMIN=EVALMIN
                      ENDIF

!                   IF (J<MAXPRINTOUT) THEN  ! COMMENTED BY BS360
                         IF (TSCONVERGED) THEN
                              WRITE(*,'(I5)',ADVANCE='NO') ITDONE
                         ELSE
                              WRITE(*,'(A5)',ADVANCE='NO') '   :('
                         ENDIF
!                   ENDIF

                    IF (ASSOCIATED(DUMMY%NEXT)) THEN
                         DUMMY=>DUMMY%NEXT
                    ELSE
                         EXIT
                    ENDIF
               ENDDO
               CALL MYCPU_TIME(ENDTIME,.FALSE.)

               WRITE(*,'(A)') '.'

               WRITE(INTSTR,'(I10)') NTSFOUND

               IF (MECCANOT) THEN                  
                  WRITE(METHSTR,'(A)') 'MECCANO'
               ELSE IF (GROWSTRINGT) THEN
                  IF (EVOLVESTRINGT) THEN
                     WRITE(METHSTR,'(A)') 'ES'
                  ELSE
                     WRITE(METHSTR,'(A)') 'GS'
                  ENDIF
               ELSE
                  WRITE(METHSTR,'(A)') 'DNEB'
               ENDIF              

               WRITE(*, '(1X,A,F7.2)',ADVANCE='YES') TRIM(METHSTR)//' RUN YIELDED '//TRIM(ADJUSTL(INTSTR))// &
                            &' TRUE TRANSITION STATE(S) TIME=',ENDTIME-STARTTIME
!              IF (NTSFOUND==1) THEN
!                   WRITE(*, '(A)') '.'
!              ELSE
!                   WRITE(*, '(A)') 'S.'
!              ENDIF
          ENDIF
          IF (SAVECANDIDATES) THEN
               IF (ASSOCIATED(FIRST)) THEN
                    DUMMY=>FIRST
                    IF (DESMINT) THEN
                       INQUIRE(IOLENGTH=RECLEN) (XYZ(3*NATOMS*(DUMMY%I-1)+1:3*NATOMS*DUMMY%I))
                    ELSE                       
                       INQUIRE(IOLENGTH=RECLEN) (XYZ(NOPT*(DUMMY%I-1)+1:NOPT*DUMMY%I))
                    ENDIF
                    J=1
                    DO
                         WRITE(FILENAME,'(I10)') J
                         FILENAME='POINTS'//TRIM(ADJUSTL(FILENAME))//'.OUT'
                         OPEN(UNIT=40,FILE=FILENAME,STATUS='UNKNOWN',FORM='UNFORMATTED',ACCESS='DIRECT',RECL=RECLEN)

                         IF (DESMINT) THEN
                            WRITE(40,REC=1) ( XYZ(3*NATOMS*(DUMMY%I-1)+1:3*NATOMS*DUMMY%I) )
                         ELSE
                            WRITE(40,REC=1) ( XYZ(NOPT*(DUMMY%I-1)+1:NOPT*DUMMY%I) )
                         ENDIF

                         CLOSE(40)
                         IF (ASSOCIATED(DUMMY%NEXT)) THEN
                              DUMMY=>DUMMY%NEXT
                              J=J+1
                         ELSE
                              EXIT
                         ENDIF
                    ENDDO
               ENDIF
          ENDIF

          IF (.NOT.ASSOCIATED(FIRST)) RETURN
          DO
               IF (.NOT.ASSOCIATED(FIRST%NEXT)) THEN
                    NULLIFY(DUMMY)
                    DEALLOCATE(FIRST)
                    RETURN
               ENDIF
               DUMMY=>FIRST%NEXT
               NULLIFY(FIRST%NEXT)
               DEALLOCATE(FIRST)
               FIRST=>DUMMY
          ENDDO
      END SUBROUTINE TSLOCATOR

SUBROUTINE CONTSLOCATOR
USE KEY,ONLY: BFGSTST,UNRST,NSTEPS,MACHINE, GROWSTRINGT, INTEPSILON, REDOTSIM
USE GSDATA, ONLY: EVOLVESTRINGT
USE KEYOUTPUT
USE MODCHARMM
USE NEBDATA
USE KEYNEB,ONLY:NIMAGE,DEBUG
USE NEBTOCONNECT
USE CHARUTILS
USE MODGUESS
USE MODMEC
USE LINKEDLIST
USE MODEFOL
USE INTCOMMONS, ONLY : DESMINT, NINTC, NNZ, KD, INTNEWT
USE COMMONS, ONLY : REDOPATH, REDOPATHNEB
IMPLICIT NONE
          
INTEGER :: I,J,ITDONE=0,J1,RECLEN,J2,MYTSMAX,NTS
INTEGER,PARAMETER :: MAXPRINTOUT = 50, ITMAX  = 30
DOUBLE PRECISION :: EDUMMY,EVALMIN,EVALMAX,MAXE,VECSNORM
LOGICAL :: TSCONVERGED
DOUBLE PRECISION,DIMENSION(3*NATOMS) :: LGDUMMY, VECS, DIAG, XLOCAL
DOUBLE PRECISION ELOCAL(NIMAGE+2)
INTEGER :: MLOC
DOUBLE PRECISION :: TIME, TIME0
DOUBLE PRECISION :: DPRAND
LOGICAL :: KNOWE, KNOWG, KNOWH ! JMC
COMMON /KNOWN/ KNOWE, KNOWG, KNOWH ! JMC
CHARACTER(LEN=256) :: FILENAME, METHSTR
LOGICAL :: TMPINTNEWT, FAILED
DOUBLE PRECISION, ALLOCATABLE :: TSGUESS(:,:), TSTEMP(:,:), LTANVEC(:,:)

MYTSMAX=10
IF (ALLOCATED(TSGUESS)) DEALLOCATE(TSGUESS)
IF (ALLOCATED(LTANVEC)) DEALLOCATE(LTANVEC)
ALLOCATE(TSGUESS(MYTSMAX,3*NATOMS),LTANVEC(MYTSMAX,3*NATOMS))
NTS=0

IF (REDOPATHNEB) THEN
   PRINT '(A,F20.10)',' CONTSLOCATOR> ERROR *** REDOPATH CANNOT BE SET WITH NEBCONSTRAINT'
   STOP
ELSE
   DO I=1,NIMAGE+1
      DO J2=1,NIMAGE+2 ! EXTRA INTERPOLATION USING THE SAME NUMBER OF IMAGES
         XLOCAL(1:NOPT)=( (NIMAGE+2-J2)*XYZ((I-1)*NOPT+1:I*NOPT)+(J2-1)*XYZ(I*NOPT+1:(I+1)*NOPT) )/(NIMAGE+1)
         CALL POTENTIAL(XLOCAL,ELOCAL(J2),LGDUMMY,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
         PRINT '(3(A,I6),A,G20.10)',' CONTSLOCATOR> ENERGY AT POSITION ',J2,' BETWEEN IMAGES ',I,' AND ',I+1, &
  &                                  ' E=',ELOCAL(J2)
      ENDDO
      IF (ELOCAL(2).LT.ELOCAL(1)) THEN
         NTS=NTS+1
         IF (NTS.GT.MYTSMAX) THEN ! INCREASE STORAGE AS REQUIRED FOR TS CANDIDATES
            ALLOCATE(TSTEMP(MYTSMAX,3*NATOMS))
            TSTEMP(1:MYTSMAX,1:3*NATOMS)=TSGUESS(1:MYTSMAX,1:3*NATOMS)
            DEALLOCATE(TSGUESS)
            ALLOCATE(TSGUESS(2*MYTSMAX,3*NATOMS))
            TSGUESS(1:MYTSMAX,1:3*NATOMS)=TSTEMP(1:MYTSMAX,1:3*NATOMS)
            TSTEMP(1:MYTSMAX,1:3*NATOMS)=LTANVEC(1:MYTSMAX,1:3*NATOMS)
            DEALLOCATE(LTANVEC)
            ALLOCATE(LTANVEC(2*MYTSMAX,3*NATOMS))
            LTANVEC(1:MYTSMAX,1:3*NATOMS)=TSTEMP(1:MYTSMAX,1:3*NATOMS)
            DEALLOCATE(TSTEMP)
            MYTSMAX=2*MYTSMAX
         ENDIF
         PRINT '(3(A,I6),A,G20.10)',' CONTSLOCATOR> ADDING TS CANDIDATE AT POSITION ',1,' BETWEEN IMAGES ',I,' AND ',I+1, &
  &                               ' E=',ELOCAL(1)
         TSGUESS(NTS,1:3*NATOMS)=XYZ((I-1)*NOPT+1:I*NOPT)
         LTANVEC(NTS,1:3*NATOMS)=XYZ((I-1)*NOPT+1:I*NOPT)-XYZ(I*NOPT+1:(I+1)*NOPT)
      ENDIF
      DO J2=2,NIMAGE+1 
         IF ( (ELOCAL(J2-1).LT.ELOCAL(J2)) .AND. (ELOCAL(J2).GT.ELOCAL(J2+1)) ) THEN
            NTS=NTS+1
            IF (NTS.GT.MYTSMAX) THEN ! INCREASE STORAGE AS REQUIRED FOR TS CANDIDATES
               ALLOCATE(TSTEMP(MYTSMAX,3*NATOMS))
               TSTEMP(1:MYTSMAX,1:3*NATOMS)=TSGUESS(1:MYTSMAX,1:3*NATOMS)
               DEALLOCATE(TSGUESS)
               ALLOCATE(TSGUESS(2*MYTSMAX,3*NATOMS))
               TSGUESS(1:MYTSMAX,1:3*NATOMS)=TSTEMP(1:MYTSMAX,1:3*NATOMS)
               TSTEMP(1:MYTSMAX,1:3*NATOMS)=LTANVEC(1:MYTSMAX,1:3*NATOMS)
               DEALLOCATE(LTANVEC)
               ALLOCATE(LTANVEC(2*MYTSMAX,3*NATOMS))
               LTANVEC(1:MYTSMAX,1:3*NATOMS)=TSTEMP(1:MYTSMAX,1:3*NATOMS)
               DEALLOCATE(TSTEMP)
               MYTSMAX=2*MYTSMAX
            ENDIF
            PRINT '(3(A,I6),A,G20.10)',' CONTSLOCATOR> ADDING TS CANDIDATE AT POSITION ',J2,' BETWEEN IMAGES ',I,' AND ',I+1, &
  &                                  ' E=',ELOCAL(J2)
            TSGUESS(NTS,1:3*NATOMS)=( (NIMAGE+2-J2)*XYZ((I-1)*NOPT+1:I*NOPT)+(J2-1)*XYZ(I*NOPT+1:(I+1)*NOPT) )/(NIMAGE+1)
            LTANVEC(NTS,1:3*NATOMS)=XYZ((I-1)*NOPT+1:I*NOPT)-XYZ(I*NOPT+1:(I+1)*NOPT)
         ENDIF
      ENDDO
      IF (ELOCAL(NIMAGE+1).LT.ELOCAL(NIMAGE+2)) THEN
         NTS=NTS+1
         IF (NTS.GT.MYTSMAX) THEN ! INCREASE STORAGE AS REQUIRED FOR TS CANDIDATES
            ALLOCATE(TSTEMP(MYTSMAX,3*NATOMS))
            TSTEMP(1:MYTSMAX,1:3*NATOMS)=TSGUESS(1:MYTSMAX,1:3*NATOMS)
            DEALLOCATE(TSGUESS)
            ALLOCATE(TSGUESS(2*MYTSMAX,3*NATOMS))
            TSGUESS(1:MYTSMAX,1:3*NATOMS)=TSTEMP(1:MYTSMAX,1:3*NATOMS)
            TSTEMP(1:MYTSMAX,1:3*NATOMS)=LTANVEC(1:MYTSMAX,1:3*NATOMS)
            DEALLOCATE(LTANVEC)
            ALLOCATE(LTANVEC(2*MYTSMAX,3*NATOMS))
            LTANVEC(1:MYTSMAX,1:3*NATOMS)=TSTEMP(1:MYTSMAX,1:3*NATOMS)
            DEALLOCATE(TSTEMP)
            MYTSMAX=2*MYTSMAX
         ENDIF
         PRINT '(3(A,I6),A,G20.10)',' CONTSLOCATOR> ADDING TS CANDIDATE AT POSITION ',NIMAGE+2,' BETWEEN IMAGES ',I,' AND ',I+1, &
  &                               ' E=',ELOCAL(NIMAGE+2)
         TSGUESS(NTS,1:3*NATOMS)=XYZ(I*NOPT+1:(I+1)*NOPT) 
         LTANVEC(NTS,1:3*NATOMS)=XYZ((I-1)*NOPT+1:I*NOPT)-XYZ(I*NOPT+1:(I+1)*NOPT)
      ENDIF
   ENDDO
ENDIF

IF (NTS.EQ.0) THEN
   PRINT '(A)',' CONTSLOCATOR> NO TS CANDIDATES TO OPTIMISE'
   STOP
ENDIF

WRITE(*,'(1X,A)',ADVANCE='NO') 'CONVERGED TO TS (NUMBER OF ITERATIONS):     '

NTSFOUND=0
CALL MYCPU_TIME(STARTTIME,.FALSE.)
DO J=1,NTS
   CALL MYCPU_TIME(TIME0,.FALSE.)
   KNOWE=.FALSE.
   KNOWG=.FALSE.
   IF (BFGSTST) THEN
      IF (UNRST) THEN 
         PRINT '(A)',' CONTSLOCATOR> ERROR *** NOT CODED FOR UNRES'
         STOP
      ELSE
         VECS(1:NOPT)=LTANVEC(J,1:NOPT)
         CALL BFGSTS(NSTEPS,TSGUESS(J,1:NOPT),EDUMMY,LGDUMMY,TSCONVERGED,RMS,EVALMIN,EVALMAX,VECS,ITDONE,.TRUE.,PRINTOPTIMIZETS)
      ENDIF
   ELSE
      CALL EFOL(TSGUESS(J,1:NOPT),TSCONVERGED,NSTEPS,EDUMMY,ITDONE,EVALMIN,DEBUG,DIAG,2)
   ENDIF
   CALL MYCPU_TIME(TIME,.FALSE.)

   IF (TSCONVERGED) THEN
      NTSFOUND=NTSFOUND+1
      ALLOCATE(TSFOUND(NTSFOUND)%E,TSFOUND(NTSFOUND)%COORD(NOPT),TSFOUND(NTSFOUND)%EVALMIN,TSFOUND(NTSFOUND)%VECS(NOPT))
      TSFOUND(NTSFOUND)%VECS=VECS(1:NOPT)
      TSFOUND(NTSFOUND)%COORD=TSGUESS(J,1:NOPT)
      TSFOUND(NTSFOUND)%E=EDUMMY
      TSFOUND(NTSFOUND)%EVALMIN=EVALMIN
      WRITE(*,'(I5)',ADVANCE='NO') ITDONE
   ELSE
      WRITE(*,'(A5)',ADVANCE='NO') '   :('
   ENDIF
ENDDO
CALL MYCPU_TIME(ENDTIME,.FALSE.)

WRITE(*,'(A)') '.'

WRITE(INTSTR,'(I10)') NTSFOUND

WRITE(*, '(A,F7.2)',ADVANCE='YES') ' CONSTRAINED POTENTIAL RUN YIELDED '//TRIM(ADJUSTL(INTSTR))// &
                  &' TRUE TRANSITION STATE(S) TIME=',ENDTIME-STARTTIME

DEALLOCATE(TSGUESS)
RETURN

END SUBROUTINE CONTSLOCATOR

      SUBROUTINE CHECKTS(DUMMY,EVALMIN,TSCONVERGED)
          USE NEBDATA
          USE MODCHARMM
          USE LINKEDLIST
          IMPLICIT NONE

          LOGICAL :: FAILCHECK,TSCONVERGED
          TYPE(CHAIN),POINTER :: DUMMY
          DOUBLE PRECISION :: EVALMIN

          ! DAE IF EVALMIN LARGE IN MAGNITUDE, THIS TS IS LIKELY TO BE BOGUS, AND CAUSE PROBLEMS
          ! WHEN THEN THE CONNECTED MINIMA HAVE TO BE FOUND
!          IF (CHRMMT.AND.(EVALMIN.LT.-100.D0)) THEN
!               TSCONVERGED=.FALSE.
!               WRITE(*,'(A,F20.10,A)') 'CHECKTS> EIGENVALUE ',EVALMIN,' TOO NEGATIVE, TS SEARCH FAILED'
!               ! DAE FOR CHARMM CHECK THIS TRANSITION STATE TO SEE IF ITS GEOMETRY HAS BECOME UNFEASIBLE
!!               CALL CHECKPOINT(XYZ(NOPT*(DUMMY%I-1):NOPT*DUMMY%I),FAILCHECK)
!!               CALL CHECKPOINT(XYZ(NOPT*(DUMMY%I-1)+1:NOPT*DUMMY%I),FAILCHECK) ! BS360
!!               IF (FAILCHECK) THEN
!!                    WRITE(*,'(A)') 'CHECKTS> TRANSITION STATE HAS UNPHYSICAL GEOMETRY, TS SEARCH FAILED'
!!                    TSCONVERGED=.FALSE.
!!               ENDIF
!          ENDIF
      END SUBROUTINE CHECKTS
END MODULE NEBOUTPUT
