!   PATHSAMPLE: A driver for OPTIM to create stationary point databases using discrete path sampling and perform kinetic analysis
!   Copyright (C) 1999-2009 David J. Wales
!   This file is part of PATHSAMPLE. 
!   
!   PATHSAMPLE is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or 
!   (at your option) any later version.
!   
!   PATHSAMPLE is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!   
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!   

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  This subroutine analyses a path.info in the new min-sad-min min-sad-min format
C  as generated with OPTIM keyword DUMPALLPATHS.
C
      SUBROUTINE GETALLPATHS
      USE PORFUNCS
      USE KEY
      USE COMMONS
      IMPLICIT NONE

      INTEGER J1, J2, ISTAT, NMINOLD, TSNUMBER, J3, NCOUNT, NMINSAVE, NTSSAVE
      DOUBLE PRECISION LOCALPOINTS(NR), ENERGY, NEWEMIN, NEWETS, DISTANCE, RMAT(3,3),
     1                 LPOINTSTS(NR), LPLUS(NR), LMINUS(NR), LOCALPOINTS2(NR)
      DOUBLE PRECISION DUMMY, DIST2 
      DOUBLE PRECISION NEWFVIBMIN, NEWFVIBTS, NEWNEGEIG, NEWPOINTSMIN(NR), NEWPOINTSMINPLUS(NR), EPLUS,
     1                 NEWPOINTSTS(NR), NEWIXMIN,  NEWIYMIN, NEWIZMIN, IXPLUS, IYPLUS, IZPLUS,
     2                 NEWIXTS,  NEWIYTS, NEWIZTS, IXMINUS, IYMINUS, IZMINUS, FRICTIONFAC
      INTEGER NEWHORDERMIN, NEWHORDERTS, NEWMIN, NEWTS, NTRIPLES
      LOGICAL TSISOLD, FAILED, MINPOLD, MINMOLD, BADTRIPLE
      CHARACTER(LEN=1) DUMMYSTRING

      NMINOLD=NMIN

      IF (MACHINE) THEN
           OPEN(1,FILE='path.info',STATUS='OLD',form='unformatted')
      ELSE
           OPEN(1,FILE='path.info',STATUS='OLD')
      ENDIF
C
C  Nasty things can happen if an OPTIM job crashes and leaves a path.info file incomplete.
C  A transition state could be written without the necessary minima. 
C  To avoid this problem, we now parse the path.info file first and only read in complete triples.
C
C  To avoid non-Morse points and higher index saddles getting into the database we now skip
C  a triple where any stationary point has an inappropriate normal mode frequency.
C  We must therefore avoid writing to the open min.data and ts.data files until we
C  have parsed the data!
C
      NCOUNT=0
      DO
         READ(1,*,END=123) DUMMYSTRING
         NCOUNT=NCOUNT+1
      ENDDO
123   REWIND(1)
      IF (NOFRQS) THEN
         NTRIPLES=NCOUNT/(3*(NATOMS+2))
         J1=NTRIPLES*3*(NATOMS+2)
         IF (DEBUG) PRINT '(2(A,I8))','getallpaths> number of triples=',NTRIPLES,' number of trailing lines=',J1-NCOUNT
      ELSE
         NTRIPLES=NCOUNT/(3*(2*NATOMS-NGLY+2))
         J1=NTRIPLES*3*(2*NATOMS-NGLY+2)
         IF (DEBUG) PRINT '(2(A,I8))','getallpaths> number of triples=',NTRIPLES,' number of trailing lines=',J1-NCOUNT
      ENDIF

      TSISOLD=.TRUE.
      DO J1=1,NTRIPLES 
         IF (DEBUG) PRINT '(A,I6,A,2I6)','getallpaths> doing triple number ',J1,' number of minima and ts=',NMIN,NTS
         BADTRIPLE=.FALSE.
C
C  NMIN and NTS can be incremented locally within the loop, but are reset to
C  NMINSAVE and NTSSAVE if we diagnose a bad triple due to bad frequencies or the
C  threshold test on the transition state fails.
C
         NMINSAVE=NMIN
         NTSSAVE=NTS
C
C  Read data for plus minimum.
C
         IF (MACHINE) THEN
            READ(1,END=110) ENERGY
            READ(1) HORDER
         ELSE
            READ(1,*,END=110) ENERGY, HORDER
         ENDIF
         NEWEMIN=ENERGY
         IF (NOFRQS) THEN
            DUMMY=4.675754133D0 ! 2 ln(2pi) +1
         ELSE
            IF (MACHINE) THEN
               READ(1) (FRQS(J2),J2=1,3*(NATOMS-NGLY))
            ELSE
               READ(1,*) (FRQS(J2),J2=1,3*(NATOMS-NGLY))
            ENDIF
            DUMMY=0.0D0
            DO J2=NFSTART,NFFINISH
               IF (FRQS(J2).LE.EVCUT) THEN
                  PRINT '(A,I8,A,G20.10)','getallpaths> SKIPPING - vibrational frequency ',J2,' of + minimum is ',FRQS(J2)
                  BADTRIPLE=.TRUE.
               ELSE
                  DUMMY=DUMMY+LOG(FRQS(J2))
               ENDIF
            ENDDO
         ENDIF
         NEWHORDERMIN=HORDER
         NEWFVIBMIN=DUMMY

         IF (MACHINE) THEN
            READ(1) (NEWPOINTSMIN(J2),J2=1,NR)  
         ELSE
            READ(1,*) (NEWPOINTSMIN(J2),J2=1,NR)  
         ENDIF
         LOCALPOINTS(1:NR)=NEWPOINTSMIN(1:NR)
         CALL INERTIAWRAPPER(LOCALPOINTS,NATOMS,angleAxis,NEWIXMIN,NEWIYMIN,NEWIZMIN)
         MINPOLD=.TRUE.
         DO J2=1,NMIN
            DISTANCE=1.0D100
            IF (ABS(NEWEMIN-EMIN(J2)).LT.EDIFFTOL) THEN
               READ(UMIN,REC=J2) (LOCALPOINTS2(J3),J3=1,NR)  
               CALL MINPERMDIST(LOCALPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY, 
     &                          RMAT,.FALSE.)
            ENDIF

            IF ((ABS(NEWEMIN-EMIN(J2)).LT.EDIFFTOL).AND.(DISTANCE.LT.GEOMDIFFTOL)) THEN
               NEWMIN=J2
               IF (DEBUG) PRINT '(2(A,I6))','getallpaths> path minimum ',2*(J1-1)+1,' is database minimum ',J2
               IF (ABS(NEWFVIBMIN-FVIBMIN(J2))/FVIBMIN(J2).GT.1.0D-3) THEN
                  WRITE(*,'(A,F15.5,A,F15.5)') 'getallpaths> WARNING, NEWFVIBMIN=',NEWFVIBMIN,' should be ',FVIBMIN(J2)
               ENDIF
               IF (NEWHORDERMIN.NE.HORDERMIN(J2)) THEN
                  WRITE(*,'(A,I6,A,I6)') 'getallpaths> ERROR, NEWHORDERMIN=',NEWHORDERMIN,' should be ',HORDERMIN(J2)
                  NEWHORDERMIN=MAX(NEWHORDERMIN,HORDERMIN(J2))
                  WRITE(*,'(A,I6)') 'getallpaths> using maximum value: ',NEWHORDERMIN
               ENDIF
               GOTO 130
            ENDIF
         ENDDO
         MINPOLD=.FALSE.
         NMIN=NMIN+1
         IF (NMIN.GT.MAXMIN) CALL MINDOUBLE
         NEWMIN=NMIN
         EMIN(NMIN)=NEWEMIN
         FVIBMIN(NMIN)=NEWFVIBMIN
         HORDERMIN(NMIN)=NEWHORDERMIN
         IF (ENSEMBLE.EQ.'T') THEN
            PFMIN(NMIN) = -EMIN(NMIN)/TEMPERATURE - FVIBMIN(NMIN)/2.0D0 - LOG(1.0D0*HORDERMIN(NMIN))
         ELSEIF (ENSEMBLE.EQ.'E') THEN
            IF (TOTALE.GT.EMIN(NMIN)) THEN
               PFMIN(NMIN) = (KAPPA-1)*LOG(TOTALE-EMIN(NMIN)) - FVIBMIN(NMIN)/2.0D0 - LOG(1.0D0*HORDERMIN(NMIN))
            ELSE
               PFMIN(NMIN) = -1.0D250
            ENDIF
         ENDIF
         PFMIN(NMIN)=PFMIN(NMIN)-PFMEAN
         IXMIN(NMIN)=NEWIXMIN
         IYMIN(NMIN)=NEWIYMIN
         IZMIN(NMIN)=NEWIZMIN
         GPFOLD(NMIN)=0.0D0
         IF (DEBUG) WRITE(*,'(A,I6,A)') 'getallpaths> new minimum ',NMIN
         TOPPOINTER(NMIN)=0
!        KSUM(NMIN)=0.0D0
C
C  We must delay this write until we know that all the stationary points are OK.
C  Writing the points to points.min doesn't matter - this is a direct access file!
C
C        WRITE(UMINDATA,'(2F20.10,I6,3F20.10)') EMIN(NMIN), FVIBMIN(NMIN), HORDERMIN(NMIN), IXMIN(NMIN), IYMIN(NMIN), IZMIN(NMIN)
C        CALL FLUSH(UMINDATA,ISTAT)
C
         WRITE(UMIN,REC=NMIN) (NEWPOINTSMIN(J2),J2=1,NR)
         CALL FLUSH(UMIN,ISTAT)

         IF (DIJINITT) THEN
            PAIRDIST(NMIN*(NMIN+1)/2)=0.0D0
            DO J3=1,NMIN-1
               READ(UMIN,REC=J3) (LOCALPOINTS(J2),J2=1,NR)
               CALL MINPERMDIST(LOCALPOINTS,NEWPOINTSMIN,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY, 
     &                          RMAT,.FALSE.)
               IF (INTERPCOSTFUNCTION) CALL MINPERMDIST(LOCALPOINTS,NEWPOINTSMIN,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,
     &                                 DISTANCE,DIST2,RIGIDBODY,RMAT,INTERPCOSTFUNCTION)
               PAIRDIST(NMIN*(NMIN-1)/2+J3)=DISTANCE
            ENDDO
         ENDIF

130      CONTINUE
         NEWPOINTSMINPLUS(1:NR)=NEWPOINTSMIN(1:NR)
         EPLUS=NEWEMIN
C
C  Read TS data.
C
         IF (MACHINE) THEN
              READ(1) ENERGY
              READ(1) HORDER
         ELSE
              READ(1,*) ENERGY, HORDER
         ENDIF
         NEWETS=ENERGY
         IF (NOFRQS) THEN
            DUMMY=1.0D0
            NEWNEGEIG=-1.0D0
         ELSE
            IF (MACHINE) THEN
               READ(1) (FRQS(J2),J2=1,3*(NATOMS-NGLY))
            ELSE
               READ(1,*) (FRQS(J2),J2=1,3*(NATOMS-NGLY))
            ENDIF
            DUMMY=0.0D0
!
! The current tests do not allow higher index saddles to get into ts.data.
!
            DO J2=NFSTART,NFFINISH-1
               IF (FRQS(J2).LT.EVCUT) THEN
                  WRITE(*,'(A,I6,A,G20.10)') 'getallpaths> SKIPPING - eigenvalue ',J2,' of this transition state is only ',FRQS(J2)
                  BADTRIPLE=.TRUE.
               ELSE
                  DUMMY=DUMMY+LOG(FRQS(J2))
               ENDIF
            ENDDO
            NEWNEGEIG=FRQS(3*(NATOMS-NGLY))
         ENDIF
C
C  Now we store the transition state coordinates.
C
         IF (MACHINE) THEN
            READ(1) (NEWPOINTSTS(J2),J2=1,NR)  
         ELSE
            READ(1,*) (NEWPOINTSTS(J2),J2=1,NR)  
         ENDIF
         LOCALPOINTS(1:NR)=NEWPOINTSTS(1:NR)
         CALL INERTIAWRAPPER(LOCALPOINTS,NATOMS,ANGLEAXIS,NEWIXTS,NEWIYTS,NEWIZTS)
         NEWFVIBTS=DUMMY
         NEWHORDERTS=HORDER
C
C  Now check for new transition states.
C
         TSISOLD=.FALSE.
!        IF ((DIJINITT.OR.DIJINITFLYT).AND.(NEWETS.GT.TSTHRESH)) THEN ! reject this ts
!
! We should also reject if both barriers are > MAXBARRIER, but we don't know the minus energy
! yet. We can allow such transition states into the database, but exclude them from
! pathway analysis.
!
         IF (NEWETS.GT.TSTHRESH) THEN ! reject this ts
            PRINT '(A,I8,A,G20.10,A)','getallpaths> Transition state ',J1,' energy ',NEWETS,' lies above the threshold'
            IF (.NOT.MINPOLD) PRINT '(A,I8,A,G20.10,A)','getallpaths> New connected plus minimum will not be added'
            BADTRIPLE=.TRUE.
            GOTO 120
         ENDIF
         DO J2=1,NTS
            DISTANCE=1.0D100
            IF (ABS(NEWETS-ETS(J2)).LT.EDIFFTOL) THEN
               READ(UTS,REC=J2) (LOCALPOINTS2(J3),J3=1,NR)
               CALL MINPERMDIST(LOCALPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY, 
     &                          RMAT,.FALSE.)
            ENDIF
!           PRINT '(A,I8,3G20.10)','J2,NEWETS,ETS(J2),DISTANCE=',J2,NEWETS,ETS(J2),DISTANCE
            IF ((ABS(NEWETS-ETS(J2)).LT.EDIFFTOL).AND.(DISTANCE.LT.GEOMDIFFTOL)) THEN
               IF (DEBUG) PRINT '(2(A,I6))','getallpaths> path ts ',J1,' is database ts ',J2
               IF (ABS(NEWFVIBTS-FVIBTS(J2))/FVIBTS(J2).GT.1.0D-3) THEN
                  WRITE(*,'(A,F15.5,A,F15.5)') 'getallpaths> WARNING, NEWFVIBTS=',NEWFVIBTS,' should be ',FVIBTS(J2)
               ENDIF
               IF (NEWHORDERTS.NE.HORDERTS(J2)) THEN
                  WRITE(*,'(A,I6,A,I6)') 'getallpaths> ERROR, NEWHORDERTS=',NEWHORDERTS,' should be ',HORDERTS(J2)
                  NEWHORDERTS=MAX(NEWHORDERTS,HORDERTS(J2))
                  WRITE(*,'(A,I6)') 'getallpaths> using maximum value: ',NEWHORDERTS
               ENDIF
               TSISOLD=.TRUE.
               TSNUMBER=J2
               GOTO 120
            ENDIF
         ENDDO
         NTS=NTS+1
         IF (NTS.GT.MAXTS) CALL TSDOUBLE
         NEWTS=NTS
         ETS(NTS)=NEWETS
         FVIBTS(NTS)=NEWFVIBTS
         HORDERTS(NTS)=NEWHORDERTS
         IXTS(NTS)=NEWIXTS
         IYTS(NTS)=NEWIYTS
         IZTS(NTS)=NEWIZTS
         NEGEIG(NTS)=NEWNEGEIG
         IF (DIJKSTRAT .OR. KSHORTESTPATHST) TSATTEMPT(NTS)=0
         IF (DEBUG) THEN
            WRITE(*,'(A,I6,A)') 'getallpaths> new intermediate ts ',NTS
         ENDIF
         PLUS(NTS)=NEWMIN
         WRITE(UTS,REC=NTS) (NEWPOINTSTS(J2),J2=1,NR)
         CALL FLUSH(UTS,ISTAT)
120      CONTINUE
C
C  Read data for minus minimum.
C
         IF (MACHINE) THEN
              READ(1) ENERGY
              READ(1) HORDER
         ELSE
              READ(1,*) ENERGY, HORDER
         ENDIF
         NEWEMIN=ENERGY
         IF (NOFRQS) THEN
            DUMMY=4.675754133D0 ! 2 ln(2pi) +1
         ELSE
            IF (MACHINE) THEN
               READ(1) (FRQS(J2),J2=1,3*(NATOMS-NGLY))
            ELSE
               READ(1,*) (FRQS(J2),J2=1,3*(NATOMS-NGLY))
            ENDIF
            DUMMY=0.0D0
            DO J2=NFSTART,NFFINISH
               IF (FRQS(J2).LE.EVCUT) THEN
                  PRINT '(A,I8,A,G20.10)','getallpaths> SKIPPING - vibrational frequency ',J2,' of - minimum is ',FRQS(J2)
                  BADTRIPLE=.TRUE.
               ELSE
                  DUMMY=DUMMY+LOG(FRQS(J2))
               ENDIF
            ENDDO
         ENDIF
         NEWHORDERMIN=HORDER
         NEWFVIBMIN=DUMMY

         IF (MACHINE) THEN
              READ(1) (NEWPOINTSMIN(J2),J2=1,NR)  
         ELSE
              READ(1,*) (NEWPOINTSMIN(J2),J2=1,NR)  
         ENDIF
         LOCALPOINTS(1:NR)=NEWPOINTSMIN(1:NR)
         CALL INERTIAWRAPPER(LOCALPOINTS,NATOMS,angleAxis,NEWIXMIN,NEWIYMIN,NEWIZMIN)
         MINMOLD=.TRUE.
         DO J2=1,NMIN
            DISTANCE=1.0D100
            IF (ABS(NEWEMIN-EMIN(J2)).LT.EDIFFTOL) THEN
               READ(UMIN,REC=J2) (LOCALPOINTS2(J3),J3=1,NR)
               CALL MINPERMDIST(LOCALPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY, 
     &                          RMAT,.FALSE.)
            ENDIF

            IF ((ABS(NEWEMIN-EMIN(J2)).LT.EDIFFTOL).AND.(DISTANCE.LT.GEOMDIFFTOL)) THEN
               NEWMIN=J2
               IF (DEBUG) PRINT '(2(A,I6))','getallpaths> path minimum ',2*J1,' is database minimum ',J2
               IF (ABS(NEWFVIBMIN-FVIBMIN(J2))/FVIBMIN(J2).GT.1.0D-3) THEN
                  WRITE(*,'(A,F15.5,A,F15.5)') 'getallpaths> WARNING, NEWFVIBMIN=',NEWFVIBMIN,' should be ',FVIBMIN(J2)
               ENDIF
               IF (NEWHORDERMIN.NE.HORDERMIN(J2)) THEN
                  WRITE(*,'(A,I6,A,I6)') 'getallpaths> ERROR, NEWHORDERMIN=',NEWHORDERMIN,' should be ',HORDERMIN(J2)
                  NEWHORDERMIN=MAX(NEWHORDERMIN,HORDERMIN(J2))
                  WRITE(*,'(A,I6)') 'getallpaths> using maximum value: ',NEWHORDERMIN
               ENDIF
               GOTO 140
            ENDIF
         ENDDO
         MINMOLD=.FALSE.
         NMIN=NMIN+1
         IF (NMIN.GT.MAXMIN) CALL MINDOUBLE
         NEWMIN=NMIN
         EMIN(NMIN)=NEWEMIN
         FVIBMIN(NMIN)=NEWFVIBMIN
         HORDERMIN(NMIN)=NEWHORDERMIN
         IF (ENSEMBLE.EQ.'T') THEN
            PFMIN(NMIN) = -EMIN(NMIN)/TEMPERATURE - FVIBMIN(NMIN)/2.0D0 - LOG(1.0D0*HORDERMIN(NMIN))
         ELSEIF (ENSEMBLE.EQ.'E') THEN
            IF (TOTALE.GT.EMIN(NMIN)) THEN
               PFMIN(NMIN) = (KAPPA-1)*LOG(TOTALE-EMIN(NMIN)) - FVIBMIN(NMIN)/2.0D0 - LOG(1.0D0*HORDERMIN(NMIN))
            ELSE
               PFMIN(NMIN) = -1.0D250
            ENDIF
         ENDIF
         PFMIN(NMIN)=PFMIN(NMIN)-PFMEAN
         IXMIN(NMIN)=NEWIXMIN
         IYMIN(NMIN)=NEWIYMIN
         IZMIN(NMIN)=NEWIZMIN
         GPFOLD(NMIN)=0.0D0
         IF (DEBUG) WRITE(*,'(A,I6,A)') 'getallpaths> new minimum ',NMIN
         TOPPOINTER(NMIN)=0
!        KSUM(NMIN)=0.0D0
C
C  Delay writing data for new minimum until transition state threshold check is passed.
C  It is OK to write points data to the direct access file.
C
C        WRITE(UMINDATA,'(2F20.10,I6,3F20.10)') EMIN(NMIN), FVIBMIN(NMIN), HORDERMIN(NMIN), IXMIN(NMIN), IYMIN(NMIN), IZMIN(NMIN)
C        CALL FLUSH(UMINDATA,ISTAT)
C
         WRITE(UMIN,REC=NMIN) (NEWPOINTSMIN(J2),J2=1,NR)
         CALL FLUSH(UMIN,ISTAT)

         IF (DIJINITT) THEN
            PAIRDIST(NMIN*(NMIN+1)/2)=0.0D0
            DO J3=1,NMIN-1
               READ(UMIN,REC=J3) (LOCALPOINTS(J2),J2=1,NR)
               CALL MINPERMDIST(LOCALPOINTS,NEWPOINTSMIN,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY,
     &                          RMAT,.FALSE.)
               IF (INTERPCOSTFUNCTION) CALL MINPERMDIST(LOCALPOINTS,NEWPOINTSMIN,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,
     &                                                  DISTANCE,DIST2,RIGIDBODY,RMAT,INTERPCOSTFUNCTION)
               PAIRDIST(NMIN*(NMIN-1)/2+J3)=DISTANCE
            ENDDO
         ENDIF

140      CONTINUE
C
C  Now we know the id of the MINUS minimum we can finish the bookkeeping for this ts.
C
         IF (TSISOLD) THEN
            IF (TSNUMBER.GT.0) THEN
               FAILED=.FALSE.
C
C  Consistency check for minima.
C
               IF (ABS(EPLUS-EMIN(PLUS(TSNUMBER))).LT.EDIFFTOL) THEN
                  IF (ABS(NEWEMIN-EMIN(MINUS(TSNUMBER))).LT.EDIFFTOL) THEN
                     CALL INERTIAWRAPPER(NEWPOINTSMINPLUS,NATOMS,angleAxis,IXPLUS,IYPLUS,IZPLUS)
                     CALL INERTIAWRAPPER(NEWPOINTSMIN,NATOMS,angleAxis,IXMINUS,IYMINUS,IZMINUS)
                  ELSE
                     PRINT '(A)','getallpaths> WARNING failed consistency check for old ts '
                     PRINT '(A,2G20.10,A,I8,A,2G20.10)','getallpaths> +/- energies=',EPLUS,NEWEMIN,' ts ',TSNUMBER,
     &                                         ' minima are: ',EMIN(PLUS(TSNUMBER)),EMIN(MINUS(TSNUMBER))
                     FAILED=.TRUE.
                  ENDIF
               ELSEIF (ABS(EPLUS-EMIN(MINUS(TSNUMBER))).LT.EDIFFTOL) THEN
                  IF (ABS(NEWEMIN-EMIN(PLUS(TSNUMBER))).LT.EDIFFTOL) THEN
                     CALL INERTIAWRAPPER(NEWPOINTSMIN,NATOMS,angleAxis,IXPLUS,IYPLUS,IZPLUS)
                     CALL INERTIAWRAPPER(NEWPOINTSMINPLUS,NATOMS,angleAxis,IXMINUS,IYMINUS,IZMINUS)
                  ELSE
                     PRINT '(A)','getallpaths> WARNING failed consistency check for old ts in getallpaths'
                     PRINT '(A,2G20.10,A,I8,A,2G20.10)','getallpaths> +/- energies=',EPLUS,NEWEMIN,' ts ',TSNUMBER,
     &                                         ' minima are: ',EMIN(PLUS(TSNUMBER)),EMIN(MINUS(TSNUMBER))
                     FAILED=.TRUE.
                  ENDIF
               ELSE
                  PRINT '(A)','getallpaths> WARNING failed consistency check for old ts'
                  PRINT '(A,2G20.10,A,I8,A,2G20.10)','getallpaths> +/- energies=',EPLUS,NEWEMIN,' ts ',TSNUMBER,
     &                                      ' minima are: ',EMIN(PLUS(TSNUMBER)),EMIN(MINUS(TSNUMBER))
                  FAILED=.TRUE.
               ENDIF
               IF ((.NOT.FAILED).AND.(.NOT.BULKT)) THEN
                  IF (ABS(IXPLUS-IXMIN(PLUS(TSNUMBER))).GT.IDIFFTOL) THEN
                     PRINT '(A,2G20.10)','getallpaths> WARNING failed consistency check for inertia: ',IXPLUS,IXMIN(PLUS(TSNUMBER))
                  ENDIF
                  IF (ABS(IYPLUS-IYMIN(PLUS(TSNUMBER))).GT.IDIFFTOL) THEN
                     PRINT '(A,2G20.10)','getallpaths> WARNING failed consistency check for inertia: ',IYPLUS,IYMIN(PLUS(TSNUMBER))
                  ENDIF
                  IF (ABS(IZPLUS-IZMIN(PLUS(TSNUMBER))).GT.IDIFFTOL) THEN
                     PRINT '(A,2G20.10)','getallpaths> WARNING failed consistency check for inertia: ',IZPLUS,IZMIN(PLUS(TSNUMBER))
                  ENDIF
                  IF (ABS(IXMINUS-IXMIN(MINUS(TSNUMBER))).GT.IDIFFTOL) THEN
                    PRINT '(A,2G20.10)','getallpaths> WARNING failed consistency check for inertia: ',IXMINUS,IXMIN(MINUS(TSNUMBER))
                  ENDIF
                  IF (ABS(IYMINUS-IYMIN(MINUS(TSNUMBER))).GT.IDIFFTOL) THEN
                    PRINT '(A,2G20.10)','getallpaths> WARNING failed consistency check for inertia: ',IYMINUS,IYMIN(MINUS(TSNUMBER))
                  ENDIF
                  IF (ABS(IZMINUS-IZMIN(MINUS(TSNUMBER))).GT.IDIFFTOL) THEN
                    PRINT '(A,2G20.10)','getallpaths> WARNING failed consistency check for inertia: ',IZMINUS,IZMIN(MINUS(TSNUMBER))
                  ENDIF
               ENDIF
            ENDIF
C
C  Old ts might link new minima. This is inconsistent and is ignored.
C  Any new minima are not written to min.data, so we should reset to the saved NMIN value,
C  just to be safe.
C
C  For small systems we see lines of control characters occasionally written to
C  min.data and ts.data. Probably due to a linux nfs problem. 
C  Try slowing down the writes somehow? 
C
            NMIN=NMINSAVE
         ELSE IF (.NOT.BADTRIPLE) THEN
            MINUS(NTS)=NEWMIN
            IF (CLOSEFILEST) OPEN(UNIT=UTSDATA,FILE='ts.data',STATUS='OLD',POSITION='APPEND')
            IF (IMFRQT) THEN
               WRITE(UTSDATA,'(2F20.10,3I10,4F20.10)') ETS(NTS),FVIBTS(NTS),HORDERTS(NTS),PLUS(NTS),MINUS(NTS),
     &                                          IXTS(NTS),IYTS(NTS),IZTS(NTS),NEGEIG(NTS)
            ELSE
               WRITE(UTSDATA,'(2F20.10,3I10,3F20.10)') ETS(NTS),FVIBTS(NTS),HORDERTS(NTS),PLUS(NTS),MINUS(NTS),
     &                                          IXTS(NTS),IYTS(NTS),IZTS(NTS)
            END IF
            CALL FLUSH(UTSDATA,ISTAT)
            IF (CLOSEFILEST) CLOSE(UNIT=UTSDATA)
            PRINT '(A,I6,A)','getallpaths> writing data for new ts to ts.data'
C           PRINT '(A,2L5,2I6)','MINPOLD,MINMOLD,NMINSAVE,NMIN=',MINPOLD,MINMOLD,NMINSAVE,NMIN
            IF (NMIN-NMINSAVE.GT.0) THEN
               PRINT '(A,I6,A)','getallpaths> writing data for ',NMIN-NMINSAVE,' new minima to min.data'
               IF (CLOSEFILEST) OPEN(UNIT=UMINDATA,FILE='min.data',STATUS='UNKNOWN',POSITION='APPEND')
               DO J2=NMINSAVE+1,NMIN
                  WRITE(UMINDATA,'(2F20.10,I6,3F20.10)') EMIN(J2),FVIBMIN(J2),HORDERMIN(J2),IXMIN(J2),IYMIN(J2),IZMIN(J2)
                  CALL FLUSH(UMINDATA,ISTAT)
               ENDDO
               IF (CLOSEFILEST) CLOSE(UNIT=UMINDATA)
            ENDIF

            IF (DIJINITT) THEN
               J2=MAX(PLUS(NTS),MINUS(NTS))
               J3=MIN(PLUS(NTS),MINUS(NTS))
               PAIRDIST(J2*(J2-1)/2+J3)=0.0D0
            ENDIF
C
C  Update ts pointers.
C
            POINTERP(NTS)=-1
            POINTERM(NTS)=-1
            IF (TOPPOINTER(PLUS(NTS)).GT.0) POINTERP(NTS)=TOPPOINTER(PLUS(NTS))
            IF (TOPPOINTER(MINUS(NTS)).GT.0) POINTERM(NTS)=TOPPOINTER(MINUS(NTS))
            TOPPOINTER(PLUS(NTS))=NTS
            TOPPOINTER(MINUS(NTS))=NTS

C           PRINT '(A,7I8)','NTS,PLUS(NTS),MINUS(NTS),TOP+,TOP-,PP,PM=',
C    &               NTS,PLUS(NTS),MINUS(NTS),TOPPOINTER(PLUS(NTS)),TOPPOINTER(MINUS(NTS)),POINTERP(NTS),POINTERM(NTS)
   
C            DO J2=NTS-1,1,-1
C               IF (PLUS(J2).EQ.PLUS(NTS)) THEN
C                  POINTERP(NTS)=J2
C                  GOTO 41
C               ELSE IF (MINUS(J2).EQ.PLUS(NTS)) THEN
C                  POINTERP(NTS)=J2
C                  GOTO 41
C               ENDIF
C            ENDDO
C41          CONTINUE
C
C            DO J2=NTS-1,1,-1
C               IF (PLUS(J2).EQ.MINUS(NTS)) THEN
C                  POINTERM(NTS)=J2
C                  GOTO 42
C               ELSE IF (MINUS(J2).EQ.MINUS(NTS)) THEN
C                  POINTERM(NTS)=J2
C                  GOTO 42
C               ENDIF
C            ENDDO
C42          CONTINUE
C
C  Calculate rates.
C
            IF (ENSEMBLE.EQ.'T') THEN
               KPLUS(NTS)=LOG(1.0D0 * HORDERMIN(PLUS(NTS))  / (2.0D0 * PI*HORDERTS(NTS))) +
     1          (FVIBMIN(PLUS(NTS))  - FVIBTS(NTS)) / 2.0D0 - (ETS(NTS) - EMIN(PLUS(NTS)))/TEMPERATURE
               IF (FRICTIONT) KPLUS(NTS)=KPLUS(NTS)+LOG(FRICTIONFAC(NEGEIG(NTS)))

               KMINUS(NTS)=LOG(1.0D0 * HORDERMIN(MINUS(NTS)) / (2.0D0 * PI*HORDERTS(NTS))) +
     1          (FVIBMIN(MINUS(NTS)) - FVIBTS(NTS)) / 2.0D0 - (ETS(NTS) - EMIN(MINUS(NTS)))/TEMPERATURE
               IF (FRICTIONT) KMINUS(NTS)=KMINUS(NTS)+LOG(FRICTIONFAC(NEGEIG(NTS)))
            ELSE
               IF (TEMPERATURE.GT.ETS(NTS)) THEN
                  KPLUS(NTS)  = LOG(1.0D0 * HORDERMIN(PLUS(NTS))  / (2*PI*HORDERTS(NTS))) +
     1               (FVIBMIN(PLUS(NTS))  - FVIBTS(NTS))/2 + (KAPPA-1)*LOG((TEMPERATURE-ETS(NTS))/(TEMPERATURE-EMIN(PLUS(NTS))))
                  KMINUS(NTS) = LOG(1.0D0 * HORDERMIN(MINUS(NTS)) / (2*PI*HORDERTS(NTS))) +
     1              (FVIBMIN(MINUS(NTS)) - FVIBTS(NTS))/2 + (KAPPA-1)*LOG((TEMPERATURE-ETS(NTS))/(TEMPERATURE-EMIN(MINUS(NTS))))
               ELSE
                  KPLUS(NTS)=-1.0D250
                  KMINUS(NTS)=-1.0D250
               ENDIF
            ENDIF
            IF (ZSYM(1:2).EQ.'CA') KPLUS(NTS)=KPLUS(NTS)+30.66356D0
            IF (ZSYM(1:2).EQ.'CA') KMINUS(NTS)=KMINUS(NTS)+30.66356D0
            IF (PLUS(NTS).EQ.MINUS(NTS)) KPLUS(NTS)=KPLUS(NTS)+LOG(2.0D0)
            IF (PLUS(NTS).EQ.MINUS(NTS)) KMINUS(NTS)=KMINUS(NTS)+LOG(2.0D0)
C
C  Don't update sum of rates out of the connected minima.
C  Assume KSUM will be calculated when needed.
C
!           IF (KSUM(PLUS(NTS)).EQ.0.0D0) THEN
!              IF (PLUS(NTS).NE.MINUS(NTS)) KSUM(PLUS(NTS))=KPLUS(NTS)
!           ELSE
!              IF (PLUS(NTS).NE.MINUS(NTS)) KSUM(PLUS(NTS))=LOG(EXP(KSUM(PLUS(NTS))-KMEAN) + EXP(KPLUS(NTS)-KMEAN)) + KMEAN
!           ENDIF
!           IF (KSUM(MINUS(NTS)).EQ.0.0D0) THEN
!              IF (PLUS(NTS).NE.MINUS(NTS)) KSUM(MINUS(NTS))=KMINUS(NTS)
!           ELSE
!              IF (PLUS(NTS).NE.MINUS(NTS)) KSUM(MINUS(NTS))=LOG(EXP(KSUM(MINUS(NTS))-KMEAN) + EXP(KMINUS(NTS)-KMEAN)) + KMEAN
!           ENDIF

            DO J2=1,NR
               LPOINTSTS(J2)=NEWPOINTSTS(J2)
               LPLUS(J2)=NEWPOINTSMINPLUS(J2)
               LMINUS(J2)=NEWPOINTSMIN(J2)
            ENDDO
            IF (ADDPT) CALL ADDPERM(LPOINTSTS,LPLUS,LMINUS) ! this is untested!
         ELSE
C
C  Old ts or bad triple encountered. Either way, resetting to saved NTS and NMIN values should be safe.
C
            NTS=NTSSAVE
            NMIN=NMINSAVE
         ENDIF
      ENDDO
110   CLOSE(1)
C
C  Find pathways out of the new intermediate minimum. We can;t do this before
C  because we need all the data for new minima and transition states to be assigned in the
C  right places.
C
      IF (CONNECTIONS.GT.2) THEN
         WRITE(*,'(A,I6,A)') 'getallpaths> checking for at least ',CONNECTIONS,' connections per minimum'
         DO J1=NMINOLD+1,NMIN
C           CALL TSSEARCH2(NEWMIN(J1),0)
            CALL TSSEARCH(J1,0)
         ENDDO
      ENDIF

      RETURN
      END
