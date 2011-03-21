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
!  Subroutine to provide candidate pairs of minima based on the minimum barrier
!  to a target minimum. Designed to remove traps.
!
SUBROUTINE GETUPAIR(NAVAIL,NUSED,MINS,MINF,SPOINTS,FPOINTS)
USE COMMON, ONLY: UMIN, NATOMS, DMIN1, DMIN2, NATTEMPT, NCPU, NMIN, PERMDIST, &
  &               NPAIRFRQ, PAIR1, PAIR2, NPAIRFRQ, NPAIRDONE, MAXPAIRS, DMINMAX, DEBUG, LOCATIONA, LOCATIONB, BULKT, &
  &               ZSYM, TWOD, DIRECTION, PLUS, MINUS, NMINA, NMINB, EMIN, NTS, ETS, EUNTRAPTHRESH, EINC, DEBUG, ANGLEAXIS, &
  &               TSTHRESH, TOPPOINTER, POINTERP, POINTERM, BOXLX, BOXLY, BOXLZ, RIGIDBODY, INTERPCOSTFUNCTION
USE PORFUNCS
IMPLICIT NONE
INTEGER NUSED, MINS, MINF, NAVAIL, PAIRSTODO, J1, J3, J4, CLOSEST(NMIN), MINVAL, OLDBASIN(NMIN), BASIN(NMIN), NBASIN, J2, NDONE
INTEGER MAXNEIGHBOURS, NP, NNEIGH, JDOING
INTEGER, ALLOCATABLE :: NEIGHBOURS(:), ITEMP(:)
DOUBLE PRECISION, ALLOCATABLE :: BLIST(:)
DOUBLE PRECISION SPOINTS(3*NATOMS), FPOINTS(3*NATOMS), BARRIER(NMIN), POINTS1(3*NATOMS), POINTS2(3*NATOMS), &
  &              DISTANCE, RMAT(3,3), DIST2, LOWESTTARG, &
  &              HIGHESTTS, DUMMY, ETHRESH
INTEGER, ALLOCATABLE :: VINT(:)
LOGICAL ISA(NMIN), ISB(NMIN), BASINT(NMIN), CHANGED, DONE(NMIN), MATCHED, OLDBASINT(NMIN)

ALLOCATE(VINT(DMINMAX))


10 CONTINUE
IF (NAVAIL.EQ.0) THEN
   PAIRSTODO=NCPU*NPAIRFRQ
   IF (NPAIRFRQ.LT.1) PAIRSTODO=NATTEMPT*NCPU ! just one set of pairs unless we run out
   IF (ALLOCATED(DMIN1)) DEALLOCATE(DMIN1,DMIN2)
   ALLOCATE(DMIN1(PAIRSTODO),DMIN2(PAIRSTODO),BLIST(PAIRSTODO))
   BLIST(1:PAIRSTODO)=-1.0D100
   CALL GETBARRIER(BARRIER,CLOSEST,LOWESTTARG)
   min: DO J1=1,NMIN
!
! Change to barrier divided by PE difference from lowest product minimum. 
! This is the PE analogue of the meric used by getfreepair.
!
      IF (EMIN(J1)-LOWESTTARG.NE.0.0D0) BARRIER(J1)=BARRIER(J1)/(EMIN(J1)-LOWESTTARG)
      IF (BARRIER(J1).LT.0.0D0) CYCLE min
      IF (CLOSEST(J1).EQ.0) CYCLE min ! sanity check !
      DO J3=1,NPAIRDONE
         IF ((PAIR1(J3).EQ.J1).AND.(PAIR2(J3).EQ.CLOSEST(J1))) CYCLE min ! do not repeat searches
      ENDDO
      NAVAIL=NAVAIL+1
      IF (NAVAIL.LE.PAIRSTODO) MINVAL=NAVAIL
      IF (PAIRSTODO.LE.NAVAIL) MINVAL=PAIRSTODO ! MIN function won;t compile here with NAG?!
      IF (MINVAL.GT.DMINMAX) THEN
         VINT(1:DMINMAX)=DMIN1(1:DMINMAX)
         DEALLOCATE(DMIN1)
         ALLOCATE(DMIN1(2*DMINMAX))
         DMIN1(1:DMINMAX)=VINT(1:DMINMAX)
         VINT(1:DMINMAX)=DMIN2(1:DMINMAX)
         DEALLOCATE(DMIN2)
         ALLOCATE(DMIN2(2*DMINMAX))
         DMIN2(1:DMINMAX)=VINT(1:DMINMAX)
         DMINMAX=2*DMINMAX
         DEALLOCATE(VINT)
         ALLOCATE(VINT(DMINMAX))
      ENDIF

      sortloop: DO J3=1,MINVAL ! sort to find the largest barriers to product
         IF (BARRIER(J1).GT.BLIST(J3)) THEN
            DO J4=MINVAL,J3+1,-1
               DMIN1(J4)=DMIN1(J4-1)
               DMIN2(J4)=DMIN2(J4-1)
               BLIST(J4)=BLIST(J4-1)
            ENDDO
            DMIN1(J3)=J1
            DMIN2(J3)=CLOSEST(J1)
            BLIST(J3)=BARRIER(J1)
            EXIT sortloop
         ENDIF
      ENDDO sortloop
      IF (DEBUG) PRINT '(3(A,I8),2(A,G20.10))','getupair> connection ',NAVAIL,' pair ',J1,' and ',CLOSEST(J1),' barrier=', &
  &                    BARRIER(J1)
   ENDDO min
   NAVAIL=MINVAL
   PRINT '(A)','getupair> saved minima pairs and barriers:'
   PRINT '(2I8,F20.10)',(DMIN1(J4),DMIN2(J4),BLIST(J4),J4=1,MINVAL)

!   GOTO 20 ! go with the highest product minimum !
!
! Now try to find a closer minimum in the product superbasin for each of the chosen
! target minima. This approach assumes that the superbasin analysis is relatively
! cheap compared with the cost of calculating many minimum distances. We look for
! closest minimum in the product superbasin that the minimum merges with that
! lies less than EUNTRAPTHRESH above. If there are no such minima then we stick
! with the highest product minimum, which is already saved in DMIN2.
! We have to repeat a lot of GETBARRIER here.
!
   ISA(1:NMIN)=.FALSE.
   ISB(1:NMIN)=.FALSE.
   DO J1=1,NMINA
      ISA(LOCATIONA(J1))=.TRUE.
   ENDDO
   DO J1=1,NMINB
      ISB(LOCATIONB(J1))=.TRUE.
   ENDDO
!
!  Find the lowest minimum in the product set and the highest transition state.
!
   LOWESTTARG=1.0D100
   IF (DIRECTION.EQ.'AB') THEN
      DO J1=1,NMINA
         IF (EMIN(LOCATIONA(J1)).LT.LOWESTTARG) LOWESTTARG=EMIN(LOCATIONA(J1))
      ENDDO
   ELSE
      DO J1=1,NMINB
         IF (EMIN(LOCATIONB(J1)).LT.LOWESTTARG) LOWESTTARG=EMIN(LOCATIONB(J1))
      ENDDO
   ENDIF
   PRINT '(A,G20.10)','getupair> lowest minimum in product set lies at ',LOWESTTARG
   HIGHESTTS=-1.0D100
   DO J1=1,NTS
      IF ((ETS(J1).GT.HIGHESTTS).AND.(ETS(J1).LT.TSTHRESH)) HIGHESTTS=ETS(J1)
   ENDDO
   PRINT '(A,G20.10)','getupair> highest transition state lies at ',HIGHESTTS
   ETHRESH=LOWESTTARG

   OLDBASIN(1:NMIN)=-1
   OLDBASINT(1:NMIN)=.FALSE.
   DONE(1:MINVAL)=.FALSE.
   NDONE=0
   DO
      BASIN(1:NMIN)=0
      NBASIN=0
      DO 
         CHANGED=.FALSE.
         DO J1=1,NTS
            IF (ETS(J1).LT.ETHRESH) THEN
               IF ((BASIN(PLUS(J1)).EQ.0).AND.(BASIN(MINUS(J1)).EQ.0)) THEN
                  CHANGED=.TRUE.
                  NBASIN=NBASIN+1
                  BASIN(PLUS(J1))=NBASIN
                  BASIN(MINUS(J1))=NBASIN
               ELSEIF (BASIN(PLUS(J1)).NE.BASIN(MINUS(J1))) THEN
                  CHANGED=.TRUE.
                  IF (BASIN(PLUS(J1)).EQ.0) THEN
                     BASIN(PLUS(J1))=BASIN(MINUS(J1))
                  ELSEIF (BASIN(MINUS(J1)).EQ.0) THEN
                     BASIN(MINUS(J1))=BASIN(PLUS(J1))
                  ELSE
!                    BASIN(PLUS(J1))=MIN(BASIN(PLUS(J1)),BASIN(MINUS(J1))) NAG can;t compile this line
                     IF (BASIN(MINUS(J1)).LT.BASIN(PLUS(J1))) BASIN(PLUS(J1))=BASIN(MINUS(J1))
                     BASIN(MINUS(J1))=BASIN(PLUS(J1))
                  ENDIF
               ENDIF
            ENDIF
!           PRINT '(A,I6,2F15.5,L5,I6)','NTS,ETS,ETHRESH,CHANGED,NBASIN=',NTS,ETS(J1),ETHRESH,CHANGED,NBASIN
         ENDDO
         IF (.NOT.CHANGED) EXIT
      ENDDO 
!
!  At this point all minima are assigned to superbasins.
!
      IF (DEBUG) PRINT '(A,I6)','superbasin analysis done NBASIN=',NBASIN
      BASINT(1:NBASIN)=.FALSE.
      IF (DIRECTION.EQ.'AB') THEN
         DO J1=1,NMINA
            IF (BASIN(LOCATIONA(J1)).GT.0) BASINT(BASIN(LOCATIONA(J1)))=.TRUE.
         ENDDO 
      ELSE
         DO J1=1,NMINB
            IF (BASIN(LOCATIONB(J1)).GT.0) BASINT(BASIN(LOCATIONB(J1)))=.TRUE.
         ENDDO 
      ENDIF
      DO J1=1,MINVAL ! cycle over target minima
         JDOING=DMIN1(J1)
         IF (BASIN(DMIN1(J1)).EQ.0) CYCLE
         IF (DONE(J1)) CYCLE
         IF (BASINT(BASIN(DMIN1(J1)))) THEN ! only search if the basins have just merged
!
!  Get connections for minimum JDOING.
!
            MAXNEIGHBOURS=10
            NNEIGH=0
            IF (ALLOCATED(NEIGHBOURS)) DEALLOCATE(NEIGHBOURS)
            ALLOCATE(NEIGHBOURS(MAXNEIGHBOURS))
            NP=TOPPOINTER(JDOING)  !  sets NP to the TS connected to minimum JDOING with the highest id
            DO WHILE (NP.GT.0)
               NNEIGH=NNEIGH+1
               IF (NNEIGH.GT.MAXNEIGHBOURS) THEN
                  ALLOCATE(ITEMP(MAXNEIGHBOURS))
                  ITEMP(1:MAXNEIGHBOURS)=NEIGHBOURS(1:MAXNEIGHBOURS)
                  DEALLOCATE(NEIGHBOURS)
                  ALLOCATE(NEIGHBOURS(2*MAXNEIGHBOURS))
                  NEIGHBOURS(1:MAXNEIGHBOURS)=ITEMP(1:MAXNEIGHBOURS)
                  DEALLOCATE(ITEMP)
                  MAXNEIGHBOURS=2*MAXNEIGHBOURS
                  IF (DEBUG) PRINT '(A,I8)','getupair> NEIGHBOURS array redimensioned size ',MAXNEIGHBOURS
               ENDIF
               MATCHED=.FALSE.
               IF (PLUS(NP).EQ.JDOING) THEN
                  MATCHED=.TRUE.
                  NEIGHBOURS(NNEIGH)=MINUS(NP)
                  NP=POINTERP(NP)
               ELSE IF (MINUS(NP).EQ.JDOING) THEN
                  MATCHED=.TRUE.
                  NEIGHBOURS(NNEIGH)=PLUS(NP)
                  NP=POINTERM(NP)
               ENDIF
               IF (.NOT.MATCHED) THEN
                  PRINT '(A,I6,A)','getupair minimum ',JDOING,' not matched - this should never happen'
                  STOP
               ENDIF
            ENDDO
            PRINT '(A,I8,A,I8,A)','getupair> minimum ',JDOING,' has ',NNEIGH,' neighbours'

            READ(UMIN,REC=DMIN1(J1)) POINTS1(1:3*NATOMS)
!
!  Which product minimum should we try to connect to? Try the closest in the
!  product superbasin before they merge within threshold EUNTRAPTHRESH.
!
            DUMMY=1.0D100
            min2: DO J2=1,NMIN
!              IF (OLDBASIN(J2).EQ.BASIN(DMIN1(J1))) THEN
               IF (OLDBASIN(J2).EQ.0) CYCLE
               IF (OLDBASINT(OLDBASIN(J2))) THEN ! this minimum was in a product basin in the last cycle
!                 PRINT '(A,7I6)','J1,DMIN1(J1),J2,BASIN(DMIN1(J1)),OLDBASIN(DMIN1(J1)),BASIN(J2),OLDBASIN(J2)=', &
!  &                               J1,DMIN1(J1),J2,BASIN(DMIN1(J1)),OLDBASIN(DMIN1(J1)),BASIN(J2),OLDBASIN(J2)
                  IF (EMIN(J2)-EMIN(DMIN1(J1)).LT.EUNTRAPTHRESH) THEN
                     DO J3=1,NNEIGH
                        IF (J2.EQ.NEIGHBOURS(J3)) THEN
                           IF (DEBUG) PRINT '(A,I6,A,I6,A)','getupair> minima ',JDOING,' and ',NEIGHBOURS(J3), &
  &                                           ' are already directly connected'
                           CYCLE min2 ! they are already connected!
                        ENDIF
                     ENDDO

                     DO J3=1,NPAIRDONE
                        IF ((PAIR1(J3).EQ.DMIN1(J1)).AND.(PAIR2(J3).EQ.J2)) CYCLE min2 ! do not repeat searches
                        IF ((PAIR1(J3).EQ.J2).AND.(PAIR2(J3).EQ.DMIN1(J1))) CYCLE min2 ! do not repeat searches
                     ENDDO
                     READ(UMIN,REC=J2) POINTS2(1:3*NATOMS)

                     CALL MINPERMDIST(POINTS1,POINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY, &
  &                                   RMAT,.FALSE.)
                     IF (INTERPCOSTFUNCTION) CALL MINPERMDIST(POINTS1,POINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD, &
  &                                                           DISTANCE,DIST2,RIGIDBODY,RMAT,INTERPCOSTFUNCTION)
                     IF (DISTANCE.LT.DUMMY) THEN
                        DUMMY=DISTANCE
                        PRINT '(3(A,I6),A,G20.10,A,F12.2)','getupair> changing partner for min ',DMIN1(J1),' from ',DMIN2(J1), &
  &                                                         ' to ',J2,' dist=', &
  &                                                          DUMMY,' ediff=',EMIN(J2)-EMIN(DMIN1(J1))
                        DMIN2(J1)=J2
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO min2
            DONE(J1)=.TRUE.
            NDONE=NDONE+1

         ENDIF
      ENDDO
      OLDBASIN(1:NMIN)=BASIN(1:NMIN)
      OLDBASINT(1:NMIN)=BASINT(1:NMIN)
      ETHRESH=ETHRESH+EINC
      IF (ETHRESH.GT.HIGHESTTS+EINC) EXIT 
      IF (NDONE.EQ.MINVAL) EXIT ! we don;t need to finish the superbasin analysis
   ENDDO

20 CONTINUE

   PRINT '(A,I8,A)','getupair> sorted list of ',NAVAIL,' pairs'
   PRINT '(2I8,F15.5)',(DMIN1(J1),DMIN2(J1),BLIST(J1),J1=1,NAVAIL)
   DEALLOCATE(BLIST)
   IF (NAVAIL.EQ.0) THEN
      PRINT '(A)','getupair> No more candidate pairs of minima in getupair - quit'
      STOP
   ENDIF
   NUSED=0
ENDIF
NUSED=NUSED+1
NAVAIL=NAVAIL-1
MINS=DMIN1(NUSED)
MINF=DMIN2(NUSED)
WRITE(*,'(5(A,I8))') 'getupair> connecting minima ',MINS,' and ',MINF, ' pairs used=',  &
  &    NUSED,' remaining=',NAVAIL,' total pairs=',NPAIRDONE
CALL FLUSH(6)
NPAIRDONE=NPAIRDONE+1
IF (NPAIRDONE.GT.MAXPAIRS) CALL PAIRDOUBLE
PAIR1(NPAIRDONE)=DMIN1(NUSED)
PAIR2(NPAIRDONE)=DMIN2(NUSED)
READ(UMIN,REC=MINS) SPOINTS(1:3*NATOMS)
READ(UMIN,REC=MINF) FPOINTS(1:3*NATOMS)

END SUBROUTINE GETUPAIR

!
!  Find the approximate lowest barrier from each minimum to a minimum from
!  the product region using a superbasin analysis.
!
SUBROUTINE GETBARRIER(BARRIER,CLOSEST,LOWESTTARG)
USE COMMON,ONLY : NMIN, NTS, ETS, EMIN, NMINA, NMINB, PLUS, MINUS, LOCATIONA, LOCATIONB, DIRECTION, EINC, DEBUG, TSTHRESH
IMPLICIT NONE
DOUBLE PRECISION HIGHESTTS, LOWESTTARG, ETHRESH, BARRIER(NMIN), DUMMY
INTEGER J1, BASIN(NMIN), CLOSEST(NMIN), NBASIN, J2
LOGICAL CHANGED, BASINT(NMIN)
LOGICAL ISA(NMIN), ISB(NMIN)

BARRIER(1:NMIN)=-1.0D0
CLOSEST(1:NMIN)=0
ISA(1:NMIN)=.FALSE.
ISB(1:NMIN)=.FALSE.
DO J1=1,NMINA
   ISA(LOCATIONA(J1))=.TRUE.
ENDDO
DO J1=1,NMINB
   ISB(LOCATIONB(J1))=.TRUE.
ENDDO
!
!  FInd the lowest minimum in the product set and the highest transition state.
!
LOWESTTARG=1.0D100
IF (DIRECTION.EQ.'AB') THEN
   DO J1=1,NMINA
      IF (EMIN(LOCATIONA(J1)).LT.LOWESTTARG) LOWESTTARG=EMIN(LOCATIONA(J1))
   ENDDO
ELSE
   DO J1=1,NMINB
      IF (EMIN(LOCATIONB(J1)).LT.LOWESTTARG) LOWESTTARG=EMIN(LOCATIONB(J1))
   ENDDO
ENDIF
PRINT '(A,G20.10)','getbarrier> lowest minimum in product set lies at ',LOWESTTARG
HIGHESTTS=-1.0D100
DO J1=1,NTS
   IF ((ETS(J1).GT.HIGHESTTS).AND.(ETS(J1).LT.TSTHRESH)) THEN
      HIGHESTTS=ETS(J1)
!     PRINT '(A,I6,G20.10)','getbarrier> J1,ETS(J1)=',J1,ETS(J1)
   ENDIF
ENDDO
PRINT '(A,G20.10)','getbarrier> highest transition state lies at ',HIGHESTTS
ETHRESH=LOWESTTARG

DO
   BASIN(1:NMIN)=0
   NBASIN=0
   DO 
      CHANGED=.FALSE.
      DO J1=1,NTS
         IF (ETS(J1).LT.ETHRESH) THEN
            IF ((BASIN(PLUS(J1)).EQ.0).AND.(BASIN(MINUS(J1)).EQ.0)) THEN
               CHANGED=.TRUE.
               NBASIN=NBASIN+1
               BASIN(PLUS(J1))=NBASIN
               BASIN(MINUS(J1))=NBASIN
            ELSEIF (BASIN(PLUS(J1)).NE.BASIN(MINUS(J1))) THEN
               CHANGED=.TRUE.
               IF (BASIN(PLUS(J1)).EQ.0) THEN
                  BASIN(PLUS(J1))=BASIN(MINUS(J1))
               ELSEIF (BASIN(MINUS(J1)).EQ.0) THEN
                  BASIN(MINUS(J1))=BASIN(PLUS(J1))
               ELSE
                  BASIN(PLUS(J1))=MIN(BASIN(PLUS(J1)),BASIN(MINUS(J1)))
                  BASIN(MINUS(J1))=BASIN(PLUS(J1))
               ENDIF
            ENDIF
         ENDIF
      ENDDO
      IF (.NOT.CHANGED) EXIT
   ENDDO 
!
!  At this point all minima are assigned to superbasins.
!
!  IF (DEBUG) PRINT '(A,I6)','superbasin analysis done, number of basins=',NBASIN
   BASINT(1:NBASIN)=.FALSE.
   IF (DIRECTION.EQ.'AB') THEN
      DO J1=1,NMINA
         IF (BASIN(LOCATIONA(J1)).GT.0) BASINT(BASIN(LOCATIONA(J1)))=.TRUE.
      ENDDO 
      DO J1=1,NMIN
         IF (BASIN(J1).EQ.0) CYCLE
         IF (BASINT(BASIN(J1)).AND.(.NOT.ISA(J1))) THEN
            IF (BARRIER(J1).LT.0.0D0) THEN
               BARRIER(J1)=ETHRESH-EMIN(J1)
!
!  Which product minimum should we try to connect to? Try the highest in the same superbasin.
!
               DUMMY=-1.0D100
               DO J2=1,NMINA
                  IF (BASIN(LOCATIONA(J2)).EQ.BASIN(J1)) THEN
                     IF (EMIN(LOCATIONA(J2)).GT.DUMMY) THEN
                        DUMMY=EMIN(LOCATIONA(J2))
                        CLOSEST(J1)=LOCATIONA(J2)
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
      ENDDO
   ELSE
      DO J1=1,NMINB
         IF (BASIN(LOCATIONB(J1)).GT.0) BASINT(BASIN(LOCATIONB(J1)))=.TRUE.
      ENDDO 
      DO J1=1,NMIN
         IF (BASIN(J1).EQ.0) CYCLE
         IF (BASINT(BASIN(J1)).AND.(.NOT.ISB(J1))) THEN
            IF (BARRIER(J1).LT.0.0D0) THEN
               BARRIER(J1)=ETHRESH-EMIN(J1)
!
!  Which product minimum should we try to connect to? Try the highest in the same superbasin.
!
               DUMMY=-1.0D100
               DO J2=1,NMINB
                  IF (BASIN(LOCATIONB(J2)).EQ.BASIN(J1)) THEN
                     IF (EMIN(LOCATIONB(J2)).GT.DUMMY) THEN
                        DUMMY=EMIN(LOCATIONB(J2))
                        CLOSEST(J1)=LOCATIONB(J2)
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF 
         ENDIF
      ENDDO
   ENDIF
   ETHRESH=ETHRESH+EINC
   IF (ETHRESH.GT.HIGHESTTS+EINC) EXIT 
ENDDO
! DO J1=1,NMIN
!    IF (BARRIER(J1).GT.0.0D0) PRINT '(A,I6,G20.10,I6)','getbarrier> J1,BARRIER,CLOSEST=',J1,BARRIER(J1),CLOSEST(J1) 
! ENDDO

END SUBROUTINE GETBARRIER

