!   NEB module is an implementation of the nudged elastic band method for performing double-ended pathway searches.
!   Copyright (C) 2003-2006 Semen A. Trygubenko and David J. Wales
!   This file is part of NEB module. NEB module is part of OPTIM.
!
!   OPTIM is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   OPTIM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
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
             PRINT '(A,F20.10)',' tslocator> transition state has energy ',EEE(REDOTSIM+1)
             ALLOCATE(FIRST)
             DUMMY=>FIRST 
             DUMMY%I=MLOC
             NULLIFY(DUMMY%NEXT)
          ELSE
             SELECT CASE(CANDIDATES) ! IDENTIFY TS CANDIDATES
             CASE('high')
                  NT=1
                  MAXE=-1.0D100
                  PRINT *,'NIMAGE=',NIMAGE
                  DO J1=2,NIMAGE+1
                     PRINT '(A,I6,2G20.10)','output> J1,EEE,MAXE=',J1,EEE(J1),MAXE
                     IF (EEE(J1).GT.MAXE) THEN
                        MLOC=J1
                        MAXE=EEE(J1)
                     ENDIF
                  ENDDO
                  ALLOCATE(FIRST)
                  DUMMY=>FIRST 
                  DUMMY%I=MLOC
                  NULLIFY(DUMMY%NEXT)
             CASE('all','maxim')
                  DO I=2,NIMAGE+1 
                       T=.TRUE.
                       IF (CANDIDATES=='maxim') then
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
               ! This should cope with highly asymmetric profiles, which otherwise looks monotonic
               ! until we try a huge number of images. 
               PRINT '(1x,a)', 'No maximum in profile - using highest image'
               ALLOCATE(FIRST)
               DUMMY=>FIRST
               IF (EEE(2).GT.EEE(NIMAGE+1)) THEN
                  DUMMY%I = 2
               ELSE
                  DUMMY%I = NIMAGE+1
               ENDIF
               NULLIFY(DUMMY%NEXT)
          ENDIF

!         write(*,'(1x,a)',advance='No') 'Following images are candidates for TS:'
!         do j=1,NTSmax
!              if (j<MaxPrintOut) write(*,'(i5)',advance='No') dummy%i-1
!              if (.not.associated(dummy%next)) exit
!              dummy=>dummy%next
!         enddo
!         write(*,'(a)') '.'

          ! ------ bs360 : more general printout (without MaxPrintOut) -------
          WRITE(*,'(1x,a,i4,a)',advance='No') 'Following ',nt,' images are candidates for TS:'
          DO J=1,NTSMAX
               WRITE(*,'(i5)',advance='No') dummy%i-1
               IF (.NOT.ASSOCIATED(DUMMY%NEXT)) EXIT
               DUMMY=>DUMMY%NEXT
               !msb50
          ENDDO
          PRINT *,' '
          ! ------ end bs360 ---------------------------
          
          IF (OPTIMIZETS) THEN
               WRITE(*,'(1x,a)',advance='No') 'Converged to TS (number of iterations):     '

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
                          KNOWG=.FALSE. ! Is this needed now that gdummy is set? DJW
                          VECS(1:NINTS)=TANVEC(1:NINTS,DUMMY%I-1)
                          VECSNORM=SUM(VECS(1:NINTS)**2)
                          IF (VECSNORM.EQ.0.0D0) THEN  ! Just in case TANVEC is somehow not set? e.g. for redopath !
                             IF (DEBUG) PRINT '(A)', ' output> setting random initial vector for eigenvector'
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
                             INTNEWT = .FALSE. ! linear transformation only
                             ! convert internal tangents to cartesians
                             CALL TRANSBACKDELTA(TANVEC(1:NOPT,DUMMY%I-1),VECS,XYZCART(3*NATOMS*(DUMMY%I-1)+1:3*NATOMS*DUMMY%I), &
                                  & NINTC,3*NATOMS,NNZ,KD,FAILED,DEBUG,INTEPSILON)                             
                             INTNEWT = TMPINTNEWT
                             VECSNORM=SUM(VECS(1:3*NATOMS)**2)
                             IF (VECSNORM.EQ.0.0D0) THEN  ! TANVEC ISN't set for GUESSPATH, MECCANO, UNMECCANO
                                IF (DEBUG) PRINT '(A)', ' output> setting random initial vector for eigenvector'
                                DO J1=1,3*NATOMS
                                   VECS(J1)=DPRAND()*2-1.0D0
                                ENDDO
                                CALL VECNORM(VECS,3*NATOMS)
                             ENDIF
                          ELSE
                             VECS(1:NOPT)=TANVEC(1:NOPT,DUMMY%I-1)
                             VECSNORM=SUM(VECS(1:NOPT)**2)
                             IF (VECSNORM.EQ.0.0D0) THEN  ! TANVEC ISN't set for GUESSPATH, MECCANO, UNMECCANO
                                IF (DEBUG) PRINT '(A)', ' output> setting random initial vector for eigenvector'
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

!                   IF (CHRMMT) CALL CHECKTS(DUMMY,EVALMIN,TSCONVERGED) ! this is now a dummy routine!

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

!                   if (j<MaxPrintOut) then  ! commented by bs360
                         IF (TSCONVERGED) THEN
                              WRITE(*,'(i5)',advance='No') itdone
                         ELSE
                              WRITE(*,'(a5)',advance='No') '   :('
                         ENDIF
!                   endif

                    IF (ASSOCIATED(DUMMY%NEXT)) THEN
                         DUMMY=>DUMMY%NEXT
                    ELSE
                         EXIT
                    ENDIF
               ENDDO
               CALL MYCPU_TIME(ENDTIME,.FALSE.)

               WRITE(*,'(a)') '.'

               WRITE(INTSTR,'(i10)') NTSfound

               IF (MECCANOT) THEN                  
                  WRITE(METHSTR,'(a)') 'MECCANO'
               ELSE IF (GROWSTRINGT) THEN
                  IF (EVOLVESTRINGT) THEN
                     WRITE(METHSTR,'(a)') 'ES'
                  ELSE
                     WRITE(METHSTR,'(a)') 'GS'
                  ENDIF
               ELSE
                  WRITE(METHSTR,'(a)') 'DNEB'
               ENDIF              

               WRITE(*, '(1x,a,f7.2)',advance='yes') trim(METHSTR)//' run yielded '//trim(adjustl(IntStr))// &
                            &' true transition state(s) time=',EndTime-StartTime
!              if (NTSfound==1) then
!                   write(*, '(a)') '.'
!              else
!                   write(*, '(a)') 's.'
!              endif
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
                         WRITE(FILENAME,'(i10)') j
                         FILENAME='points'//trim(adjustl(filename))//'.out'
                         OPEN(UNIT=40,FILE=FILENAME,STATUS='unknown',form='unformatted',access='direct',recl=reclen)

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
   PRINT '(A,F20.10)',' contslocator> ERROR *** REDOPATH cannot be set with NEBCONSTRAINT'
   STOP
ELSE
   DO I=1,NIMAGE+1
      DO J2=1,NIMAGE+2 ! extra interpolation using the same number of images
         XLOCAL(1:NOPT)=( (NIMAGE+2-J2)*XYZ((I-1)*NOPT+1:I*NOPT)+(J2-1)*XYZ(I*NOPT+1:(I+1)*NOPT) )/(NIMAGE+1)
         CALL POTENTIAL(XLOCAL,ELOCAL(J2),LGDUMMY,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
         PRINT '(3(A,I6),A,G20.10)',' contslocator> energy at position ',J2,' between images ',I,' and ',I+1, &
  &                                  ' E=',ELOCAL(J2)
      ENDDO
      IF (ELOCAL(2).LT.ELOCAL(1)) THEN
         NTS=NTS+1
         IF (NTS.GT.MYTSMAX) THEN ! increase storage as required for TS candidates
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
         PRINT '(3(A,I6),A,G20.10)',' contslocator> adding ts candidate at position ',1,' between images ',I,' and ',I+1, &
  &                               ' E=',ELOCAL(1)
         TSGUESS(NTS,1:3*NATOMS)=XYZ((I-1)*NOPT+1:I*NOPT)
         LTANVEC(NTS,1:3*NATOMS)=XYZ((I-1)*NOPT+1:I*NOPT)-XYZ(I*NOPT+1:(I+1)*NOPT)
      ENDIF
      DO J2=2,NIMAGE+1 
         IF ( (ELOCAL(J2-1).LT.ELOCAL(J2)) .AND. (ELOCAL(J2).GT.ELOCAL(J2+1)) ) THEN
            NTS=NTS+1
            IF (NTS.GT.MYTSMAX) THEN ! increase storage as required for TS candidates
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
            PRINT '(3(A,I6),A,G20.10)',' contslocator> adding ts candidate at position ',J2,' between images ',I,' and ',I+1, &
  &                                  ' E=',ELOCAL(J2)
            TSGUESS(NTS,1:3*NATOMS)=( (NIMAGE+2-J2)*XYZ((I-1)*NOPT+1:I*NOPT)+(J2-1)*XYZ(I*NOPT+1:(I+1)*NOPT) )/(NIMAGE+1)
            LTANVEC(NTS,1:3*NATOMS)=XYZ((I-1)*NOPT+1:I*NOPT)-XYZ(I*NOPT+1:(I+1)*NOPT)
         ENDIF
      ENDDO
      IF (ELOCAL(NIMAGE+1).LT.ELOCAL(NIMAGE+2)) THEN
         NTS=NTS+1
         IF (NTS.GT.MYTSMAX) THEN ! increase storage as required for TS candidates
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
         PRINT '(3(A,I6),A,G20.10)',' contslocator> adding ts candidate at position ',NIMAGE+2,' between images ',I,' and ',I+1, &
  &                               ' E=',ELOCAL(NIMAGE+2)
         TSGUESS(NTS,1:3*NATOMS)=XYZ(I*NOPT+1:(I+1)*NOPT) 
         LTANVEC(NTS,1:3*NATOMS)=XYZ((I-1)*NOPT+1:I*NOPT)-XYZ(I*NOPT+1:(I+1)*NOPT)
      ENDIF
   ENDDO
ENDIF

IF (NTS.EQ.0) THEN
   PRINT '(A)',' contslocator> No ts candidates to optimise'
   STOP
ENDIF

WRITE(*,'(1x,a)',advance='No') 'Converged to TS (number of iterations):     '

NTSFOUND=0
CALL MYCPU_TIME(STARTTIME,.FALSE.)
DO J=1,NTS
   CALL MYCPU_TIME(TIME0,.FALSE.)
   KNOWE=.FALSE.
   KNOWG=.FALSE.
   IF (BFGSTST) THEN
      IF (UNRST) THEN 
         PRINT '(A)',' contslocator> ERROR *** not coded for UNRES'
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
      WRITE(*,'(i5)',advance='No') itdone
   ELSE
      WRITE(*,'(a5)',advance='No') '   :('
   ENDIF
ENDDO
CALL MYCPU_TIME(ENDTIME,.FALSE.)

WRITE(*,'(a)') '.'

WRITE(INTSTR,'(i10)') NTSfound

WRITE(*, '(A,F7.2)',advance='yes') ' Constrained potential run yielded '//trim(adjustl(IntStr))// &
                  &' true transition state(s) time=',EndTime-StartTime

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

          ! DAE If EVALMIN large in magnitude, this TS is likely to be bogus, and cause problems
          ! when then the connected minima have to be found
!          IF (CHRMMT.AND.(EVALMIN.LT.-100.D0)) THEN
!               TSConverged=.FALSE.
!               WRITE(*,'(A,F20.10,A)') 'checkts> Eigenvalue ',EVALMIN,' too negative, TS search failed'
!               ! DAE for CHARMM check this transition state to see if its geometry has become unfeasible
!!               CALL CHECKPOINT(xyz(nopt*(dummy%i-1):nopt*dummy%i),FAILCHECK)
!!               CALL CHECKPOINT(xyz(nopt*(dummy%i-1)+1:nopt*dummy%i),FAILCHECK) ! bs360
!!               IF (FAILCHECK) THEN
!!                    WRITE(*,'(A)') 'checkts> Transition state has unphysical geometry, TS search failed'
!!                    TSConverged=.FALSE.
!!               ENDIF
!          ENDIF
      END SUBROUTINE CHECKTS
END MODULE NEBOUTPUT
