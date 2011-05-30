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

!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  This subroutine merges the database information in directory PATHNAME.
!
SUBROUTINE MERGEDB
USE COMMON,ONLY : NATOMS, IYTS, IZTS, UTSDATA, UTS, HORDERMIN, TOPPOINTER, HORDERTS, PLUS, MINUS, GPFOLD, &
   &              MAXMIN, MAXTS, FVIBTS, EMIN, FVIBMIN, IXMIN, IYMIN, IZMIN, NEGEIG, PATHNAME, ETS, DEBUG, &
   &              NMIN, UNRST, CHARMMT, IDIFFTOL, EDIFFTOL, UMINDATA, UMIN, NTS, IXTS, NMINA, NMINB, &
   &              LOCATIONA, LOCATIONB, ANGLEAXIS, PERMDIST, BOXLX, BOXLY, BOXLZ, GEOMDIFFTOL, TWOD, &
   &              RIGIDBODY, BULKT, ZSYM, PERMISOMER, IMFRQT, CLOSEFILEST, AMHT

USE PORFUNCS
IMPLICIT NONE

INTEGER J1, J2, ISTAT, NMINOLD, NTSOLD, NMINDB, NDUMMY, J3
DOUBLE PRECISION LOCALPOINTS(3*NATOMS), NEWEMIN, NEWETS
DOUBLE PRECISION NEWFVIBMIN, NEWFVIBTS, NEWPOINTSMIN(3*NATOMS), NEWNEGEIG, &
  &  NEWPOINTSTS(3*NATOMS), NEWIXMIN,  NEWIYMIN, NEWIZMIN, &
  &  NEWIXTS,  NEWIYTS, NEWIZTS, DISTANCE, DIST2, LOCALPOINTS2(3*NATOMS), RMAT(3,3)
INTEGER NEWHORDERMIN, NEWHORDERTS, NEWMIN, NEWTS, INDEX
INTEGER, ALLOCATABLE :: MINMAP(:), IDUM(:), TSMAP(:), LOCATIONDB(:)
CHARACTER(LEN=130) FNAME
INTEGER MINDB, TSDB
INTEGER :: NMINDBMAX=10
INTEGER :: NTSDBMAX=10

!
! First check for new minima in <pathname>/min.data and <pathname>/points.min
!
NMINOLD=NMIN
FNAME=TRIM(ADJUSTL(PATHNAME)) // '/min.data'
PRINT *,'FNAME=',FNAME
OPEN(UNIT=1,FILE=TRIM(ADJUSTL(FNAME)),STATUS='OLD')
FNAME=TRIM(ADJUSTL(PATHNAME)) // '/points.min'
OPEN(UNIT=2,FILE=TRIM(ADJUSTL(FNAME)),ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='OLD',RECL=8*3*NATOMS)
NEWMIN=0
MINDB=0
ALLOCATE(MINMAP(NMINDBMAX))
DO
   MINDB=MINDB+1
   IF (MINDB.GT.NMINDBMAX) THEN
      ALLOCATE(IDUM(NMINDBMAX))
      IDUM(1:NMINDBMAX)=MINMAP(1:NMINDBMAX)
      DEALLOCATE(MINMAP)
      ALLOCATE(MINMAP(2*NMINDBMAX))
      MINMAP(1:NMINDBMAX)=IDUM(1:NMINDBMAX)
      NMINDBMAX=2*NMINDBMAX
      DEALLOCATE(IDUM)
   ENDIF
   INDEX=NMIN+NEWMIN+1
   IF (INDEX.GT.MAXMIN) CALL MINDOUBLE
   IF (IMFRQT) THEN
      READ(1,*,END=30) EMIN(INDEX),FVIBMIN(INDEX),HORDERMIN(INDEX),IXMIN(INDEX),IYMIN(INDEX),IZMIN(INDEX),NEGEIG(INDEX)
   ELSE
      READ(1,*,END=30) EMIN(INDEX),FVIBMIN(INDEX),HORDERMIN(INDEX),IXMIN(INDEX),IYMIN(INDEX),IZMIN(INDEX)
   END IF
   NEWEMIN=EMIN(INDEX)
   NEWFVIBMIN=FVIBMIN(INDEX)
   NEWHORDERMIN=HORDERMIN(INDEX)
!
!  Read in points and check for agreement with moments of inertia as in setup
!
   READ(2,REC=MINDB) (NEWPOINTSMIN(J2),J2=1,3*NATOMS)  
   LOCALPOINTS(1:3*NATOMS)=NEWPOINTSMIN(1:3*NATOMS)
   IF (AMHT) THEN
      WRITE(*,*)'mergedb> AVOIDING MOMENT OF INITERIA CALC FOR AMH'
   ELSE
      CALL INERTIAWRAPPER(LOCALPOINTS,NATOMS,angleAxis,NEWIXMIN,NEWIYMIN,NEWIZMIN)
      IF ((ABS(NEWIXMIN-IXMIN(INDEX)).GT.IDIFFTOL).OR. &
  &    (ABS(NEWIYMIN-IYMIN(INDEX)).GT.IDIFFTOL).OR. &
  &    (ABS(NEWIZMIN-IZMIN(INDEX)).GT.IDIFFTOL)) THEN
         WRITE(*,'(A)') 'mergedb> possible error - principal moments of inertia do not agree with input'
         WRITE(*,'(A,3F20.10)') 'mergedb> values from coordinates: ',NEWIXMIN,NEWIYMIN,NEWIZMIN
         WRITE(*,'(A,3F20.10)') 'mergedb> values from ' // TRIM(ADJUSTL(PATHNAME)) // '/min.data', &
   &                                     IXMIN(INDEX),IYMIN(INDEX),IZMIN(INDEX)
        STOP
      ENDIF
   ENDIF
!
!  Is it a new minimum? Set up MINMAP.
!  Testing the new minima against themselves allows us to remove duplicates from 
!  the database we are merging! 
!
!  DO J2=1,NMIN
   DO J2=1,NMIN+NEWMIN
      IF (ABS(NEWEMIN-EMIN(J2)).LT.EDIFFTOL) THEN
         DISTANCE=1.0D100
         READ(UMIN,REC=J2) (LOCALPOINTS2(J3),J3=1,3*NATOMS)
         CALL MINPERMDIST(LOCALPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY, &
  &                       RMAT,.FALSE.)
         PRINT '(A,2F20.10)','DISTANCE,GEOMDIFFTOL=',DISTANCE,GEOMDIFFTOL
         IF ((ABS(NEWEMIN-EMIN(J2)).LT.1.0D-9).AND.(DISTANCE.GT.GEOMDIFFTOL)) THEN
            PRINT '(A)',' likely error?'
            PRINT '(A,2G20.10)',' NEWEMIN,EMIN(J2)=',NEWEMIN,EMIN(J2)
            PRINT '(A,2I6)',' J2,J3=',J2,J3
            PRINT '(A)','LOCALPOINTS:'
            PRINT '(3F20.10)',LOCALPOINTS(1:3*NATOMS)
            PRINT '(A)','LOCALPOINTS2:'
            PRINT '(3F20.10)',LOCALPOINTS2(1:3*NATOMS)
            STOP !!! DJW
         ENDIF
         
         IF (DISTANCE.LT.GEOMDIFFTOL) THEN
            IF (DEBUG) PRINT '(2(A,I6))','mergedb> minimum ',MINDB,' is database minimum ',J2
            IF (ABS(NEWFVIBMIN-FVIBMIN(J2))/FVIBMIN(J2).GT.1.0D-3) THEN
               WRITE(*,'(A,F15.5,A,F15.5)') 'mergedb> WARNING, NEWFVIBMIN=',NEWFVIBMIN,' should be ',FVIBMIN(J2)
            ENDIF
            IF (NEWHORDERMIN.NE.HORDERMIN(J2)) THEN
               WRITE(*,'(A,I6,A,I6)') 'mergedb> ERROR, NEWHORDERMIN=',NEWHORDERMIN,' should be ',HORDERMIN(J2)
               NEWHORDERMIN=HORDERMIN(J2)
               WRITE(*,'(A,I6)') 'mergedb> using existing value: ',HORDERMIN(J2)
            ENDIF
            MINMAP(MINDB)=J2
            GOTO 130
         ENDIF
      ENDIF
   ENDDO
!
!  If we reach this point we have a new minimum
!
   NEWMIN=NEWMIN+1
   MINMAP(MINDB)=NMIN+NEWMIN
   GPFOLD(NMIN+NEWMIN)=0.0D0
   TOPPOINTER(NMIN+NEWMIN)=0
   IF (DEBUG) PRINT '(2(A,I6))','mergedb> new minimum number ',NMIN+NEWMIN,' number of new minima=',NEWMIN
   IF (DEBUG) WRITE(*,'(A,I6,A)') 'mergedb> new minimum ',NMIN+NEWMIN,&
  &                ' writing parameters to file min.data and points to points.min'
   IF (CLOSEFILEST) OPEN(UNIT=UMINDATA,FILE='min.data',STATUS='UNKNOWN',POSITION='APPEND')
   WRITE(UMINDATA,'(2F20.10,I6,3F20.10)') EMIN(INDEX),FVIBMIN(INDEX),HORDERMIN(INDEX),IXMIN(INDEX),IYMIN(INDEX),IZMIN(INDEX)
   CALL FLUSH(UMINDATA,ISTAT)
   IF (CLOSEFILEST) CLOSE(UNIT=UMINDATA)
   WRITE(UMIN,REC=INDEX) (NEWPOINTSMIN(J2),J2=1,3*NATOMS)
   
130 CONTINUE
ENDDO
30 NMIN=NMINOLD+NEWMIN
CLOSE(1)
CLOSE(2)
IF (NMIN.GT.MAXMIN) CALL MINDOUBLE
MINDB=MINDB-1
PRINT '(A,I6,A,I6,A)','mergedb> merged ',MINDB,' minima - ',NEWMIN,' are new'
!
!  Now for the transition states.
!
NTSOLD=NTS
FNAME=TRIM(ADJUSTL(PATHNAME)) // '/ts.data'
OPEN(UNIT=1,FILE=TRIM(ADJUSTL(FNAME)),STATUS='OLD')
FNAME=TRIM(ADJUSTL(PATHNAME)) // '/points.ts'
OPEN(UNIT=2,FILE=TRIM(ADJUSTL(FNAME)),ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='OLD',RECL=8*3*NATOMS)
NEWTS=0
TSDB=0
NDUMMY=0
ALLOCATE(TSMAP(NTSDBMAX))
DO
   TSDB=TSDB+1
   IF (TSDB.GT.NTSDBMAX) THEN
      ALLOCATE(IDUM(NTSDBMAX))
      IDUM(1:NTSDBMAX)=TSMAP(1:NTSDBMAX)
      DEALLOCATE(TSMAP)
      ALLOCATE(TSMAP(2*NTSDBMAX))
      TSMAP(1:NTSDBMAX)=IDUM(1:NTSDBMAX)
      NTSDBMAX=2*NTSDBMAX
      DEALLOCATE(IDUM)
   ENDIF
   INDEX=NTS+NEWTS+1
   IF (INDEX.GT.MAXTS) CALL TSDOUBLE
   IF (IMFRQT) THEN
      READ(1,*,END=40) &
         ETS(INDEX),FVIBTS(INDEX),HORDERTS(INDEX),PLUS(INDEX),MINUS(INDEX),IXTS(INDEX),IYTS(INDEX),IZTS(INDEX),NEGEIG(INDEX)
   ELSE
      READ(1,*,END=40) ETS(INDEX),FVIBTS(INDEX),HORDERTS(INDEX),PLUS(INDEX),MINUS(INDEX),IXTS(INDEX),IYTS(INDEX),IZTS(INDEX)
   ENDIF
   NEWETS=ETS(INDEX)
   NEWFVIBTS=FVIBTS(INDEX)
   NEWHORDERTS=HORDERTS(INDEX)
   NEWNEGEIG=NEGEIG(INDEX)
!
!  Read in points and check for agreement with moments of inertia as in setup
!
   READ(2,REC=TSDB) (NEWPOINTSTS(J2),J2=1,3*NATOMS)  
   LOCALPOINTS(1:3*NATOMS)=NEWPOINTSTS(1:3*NATOMS)
   IF (AMHT) THEN
      WRITE(*,*)'mergedb> AVOIDING MOMENT OF INITERIA CALC FOR AMH'
   ELSE
      CALL INERTIAWRAPPER(LOCALPOINTS,NATOMS,angleAxis,NEWIXTS,NEWIYTS,NEWIZTS)
      IF ((ABS(NEWIXTS-IXTS(INDEX)).GT.IDIFFTOL).OR. &
  &       (ABS(NEWIYTS-IYTS(INDEX)).GT.IDIFFTOL).OR. &
  &       (ABS(NEWIZTS-IZTS(INDEX)).GT.IDIFFTOL)) THEN
         WRITE(*,'(A)') 'mergedb> possible error - principal moments of inertia do not agree with input'
         WRITE(*,'(A,3F20.10)') 'mergedb> values from coordinates: ',NEWIXTS,NEWIYTS,NEWIZTS
         WRITE(*,'(A,3F20.10)') 'mergedb> values from ' // TRIM(ADJUSTL(PATHNAME)) // '/ts.data', &
   &                                        IXMIN(INDEX),IYMIN(INDEX),IZMIN(INDEX)
         IF (.NOT.BULKT) STOP
      ENDIF
   ENDIF
!
!  Is it a new ts? Set up TSMAP.
!
!  DO J2=1,NTS
   DO J2=1,NTS+NEWTS
      IF (ABS(NEWETS-ETS(J2)).LT.EDIFFTOL) THEN
         DISTANCE=1.0D100
         READ(UTS,REC=J2) (LOCALPOINTS2(J3),J3=1,3*NATOMS)
         CALL MINPERMDIST(LOCALPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY, &
  &                       RMAT,.FALSE.)
         IF (DISTANCE.LT.GEOMDIFFTOL) THEN

            IF (DEBUG) PRINT '(2(A,I6))','mergedb> ts ',TSDB,' is database ts ',J2
            IF (ABS(NEWFVIBTS-FVIBTS(J2))/FVIBTS(J2).GT.1.0D-3) THEN
               WRITE(*,'(A,F15.5,A,F15.5)') 'mergedb> WARNING, NEWFVIBTS=',NEWFVIBTS,' should be ',FVIBTS(J2)
            ENDIF
            IF (NEWHORDERTS.NE.HORDERTS(J2)) THEN
            WRITE(*,'(A,I6,A,I6)') 'mergedb> ERROR, NEWHORDERTS=',NEWHORDERTS,' should be ',HORDERTS(J2)
                  NEWHORDERTS=HORDERTS(J2)
               WRITE(*,'(A,I6)') 'mergedb> using existing value: ',HORDERTS(J2)
            ENDIF
!
!  Check end minima are consistent as well.
!
            IF (.NOT.(((MINMAP(PLUS(INDEX)).EQ.PLUS(J2)).AND.(MINMAP(MINUS(INDEX)).EQ.MINUS(J2))).OR. &
  &                   ((MINMAP(PLUS(INDEX)).EQ.MINUS(J2)).AND.(MINMAP(MINUS(INDEX)).EQ.PLUS(J2)))) ) THEN
               PRINT '(A,I6)', 'mergedb> Inconsistent minima for old transition state ',J2
               PRINT '(A,2I6)','mergedb> Previous minima: ',PLUS(J2),MINUS(J2)
               PRINT '(A,2I6)','mergedb> Minima for database '//TRIM(ADJUSTL(PATHNAME))//'/ts.data: ', &
  &                             MINMAP(PLUS(INDEX)),MINMAP(MINUS(INDEX))
            ENDIF
            TSMAP(TSDB)=J2
            NDUMMY=NDUMMY+1
            GOTO 140
         ENDIF
      ENDIF
   ENDDO
!
!  If we reach this point we have a new ts
!
   NEWTS=NEWTS+1
   TSMAP(TSDB)=NTS+NEWTS
   IF (DEBUG) PRINT '(A,I6,A,I6)','mergedb> new ts number ',NTS+NEWTS,' number of new ts=',NEWTS
   IF (DEBUG) WRITE(*,'(A,I6,A)') 'mergedb> new ts ',NEWTS,' writing parameters to file ts.data and points to points.ts'
   IF (CLOSEFILEST) OPEN(UNIT=UTSDATA,FILE='ts.data',STATUS='OLD',POSITION='APPEND')
   IF (IMFRQT) THEN
      WRITE(UTSDATA,'(2F20.10,3I10,4F20.10)') ETS(INDEX),FVIBTS(INDEX),HORDERTS(INDEX),MINMAP(PLUS(INDEX)),MINMAP(MINUS(INDEX)), &
        &                                          IXTS(INDEX),IYTS(INDEX),IZTS(INDEX),NEGEIG(INDEX)
   ELSE
      WRITE(UTSDATA,'(2F20.10,3I10,3F20.10)') ETS(INDEX),FVIBTS(INDEX),HORDERTS(INDEX),MINMAP(PLUS(INDEX)),MINMAP(MINUS(INDEX)), &
        &                                          IXTS(INDEX),IYTS(INDEX),IZTS(INDEX)
   ENDIF
   CALL FLUSH(UTSDATA,ISTAT)
   IF (CLOSEFILEST) CLOSE(UNIT=UTSDATA)
   WRITE(UTS,REC=INDEX) (NEWPOINTSTS(J2),J2=1,3*NATOMS)
   
140 CONTINUE
ENDDO
40 NTS=NTSOLD+NEWTS
CLOSE(1)
CLOSE(2)
IF (NTS.GT.MAXTS) CALL TSDOUBLE
TSDB=TSDB-1
PRINT '(A,I6,A,I6,A)','mergedb> merged ',TSDB,' ts - ',NEWTS,' are new'
!
!  Finally, check min.A and min.B
!

FNAME=TRIM(ADJUSTL(PATHNAME)) // '/min.A'
OPEN(UNIT=1,FILE=TRIM(ADJUSTL(FNAME)),STATUS='OLD')
READ(1,*) NMINDB
ALLOCATE(LOCATIONDB(NMINDB)) 
READ(1,*) LOCATIONDB(1:NMINDB) 
CLOSE(1)
NEWMIN=0
OPEN(UNIT=3,FILE='min.A.new',STATUS='UNKNOWN')
DO J1=1,NMINDB
   DO J2=1,NMINA
      IF (MINMAP(LOCATIONDB(J1)).EQ.LOCATIONA(J1)) GOTO 55 ! it was already a A minimum
   ENDDO
   NEWMIN=NEWMIN+1
55 CONTINUE
ENDDO
PRINT '(A,I6)','mergedb> number of new A minima: ',NEWMIN
WRITE(3,*) NMINA+NEWMIN
IF (NMINA.GT.0) WRITE(3,'(I6)') LOCATIONA(1:NMINA)
DO J1=1,NMINDB
   PRINT '(3(A,I6))','mergedb> A minimum ',J1,' was number ',LOCATIONDB(J1),' maps to number ',MINMAP(LOCATIONDB(J1))
   DO J2=1,NMINA
      IF (MINMAP(LOCATIONDB(J1)).EQ.LOCATIONA(J1)) GOTO 50 ! it was already an A minimum
   ENDDO
   WRITE(3,'(I6)') MINMAP(LOCATIONDB(J1))
50 CONTINUE
ENDDO
CLOSE(3)
DEALLOCATE(LOCATIONDB)

FNAME=TRIM(ADJUSTL(PATHNAME)) // '/min.B'
OPEN(UNIT=1,FILE=TRIM(ADJUSTL(FNAME)),STATUS='OLD')
READ(1,*) NMINDB
ALLOCATE(LOCATIONDB(NMINDB)) 
READ(1,*) LOCATIONDB(1:NMINDB) 
CLOSE(1)
NEWMIN=0
OPEN(UNIT=3,FILE='min.B.new',STATUS='UNKNOWN')
DO J1=1,NMINDB
   DO J2=1,NMINB
      IF (MINMAP(LOCATIONDB(J1)).EQ.LOCATIONB(J1)) GOTO 65 ! it was already a B minimum
   ENDDO
   NEWMIN=NEWMIN+1
65 CONTINUE
ENDDO
WRITE(3,*) NMINB+NEWMIN
PRINT '(A,I6)','mergedb> number of new B minima: ',NEWMIN
IF (NMINB.GT.0) WRITE(3,'(I6)') LOCATIONB(1:NMINB)
DO J1=1,NMINDB
   PRINT '(3(A,I6))','mergedb> B minimum ',J1,' was number ',LOCATIONDB(J1),' maps to number ',MINMAP(LOCATIONDB(J1))
   DO J2=1,NMINB
      IF (MINMAP(LOCATIONDB(J1)).EQ.LOCATIONB(J1)) GOTO 60 ! it was already a B minimum
   ENDDO
   WRITE(3,'(I6)') MINMAP(LOCATIONDB(J1))
60 CONTINUE
ENDDO
CLOSE(3)
DEALLOCATE(LOCATIONDB)

STOP
END
