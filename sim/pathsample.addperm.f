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
C  Check transition state NTS and the minima that it connects to see if
C  permutational isomers can be added to those currently known.
C
C  One potential problem: the orders of the point groups are not
C  calculated for new stationary points. Of course, they are all one for
C  this system. This affects the rates as well, so we can.t just copy
C  the corresponding quantities from the original transition state.
C
      SUBROUTINE ADDPERM(LPOINTSTS,LPLUS,LMINUS)
      USE PORFUNCS
      USE COMMON
      IMPLICIT NONE
      INTEGER J1,J2,DOTS,ISTAT
      DOUBLE PRECISION LPOINTS(3*NATOMS), LDUMMY,NIX,NIY,NIZ,LPOINTSTS(3*NATOMS),LPLUS(3*NATOMS),
     1                 LMINUS(3*NATOMS), FRICTIONFAC

      PRINT '(A)','WARNING- routine ADDPERM has not been tested'
      IF (NTAG.LE.0) THEN
         WRITE(*,'(A,I5)') 'WARNING - addperm called but NTAG=',NTAG
         RETURN
      ENDIF
      DOTS=NTS

      DO J1=2,NATOMS
         DO J2=1,3*NATOMS
            LPOINTS(J2)=LPOINTSTS(J2)
         ENDDO
C
C  Permute tagged atom NTAG and atom J1
C
         LDUMMY=LPOINTS(3*(NTAG-1)+1)
         LPOINTS(3*(NTAG-1)+1)=LPOINTS(3*(J1-1)+1)
         LPOINTS(3*(J1-1)+1)=LDUMMY
         LDUMMY=LPOINTS(3*(NTAG-1)+2)
         LPOINTS(3*(NTAG-1)+2)=LPOINTS(3*(J1-1)+2)
         LPOINTS(3*(J1-1)+2)=LDUMMY
         LDUMMY=LPOINTS(3*(NTAG-1)+3)
         LPOINTS(3*(NTAG-1)+3)=LPOINTS(3*(J1-1)+3)
         LPOINTS(3*(J1-1)+3)=LDUMMY

         CALL INERTIAWRAPPER(LPOINTS,NATOMS,angleAxis,NIX,NIY,NIZ)
         DO J2=1,NTS
            IF ((ABS(ETS(J2)-ETS(DOTS)).LT.EDIFFTOL).AND.
     1          (ABS(IXTS(J2)-NIX).LT.IDIFFTOL).AND.
     2          (ABS(IYTS(J2)-NIY).LT.IDIFFTOL).AND.
     3          (ABS(IZTS(J2)-NIZ).LT.IDIFFTOL)) THEN
               WRITE(*,'(A,I5,A,I5,A,I5)') 
     1                 'permuting atoms 1 and ',J1,' for ts ',DOTS,' gives old transition state ',J2
               GOTO 10
            ENDIF
         ENDDO
         WRITE(*,'(A,I5,A,I5,A,I5)') 'permuting atoms 1 and ',J1,' for ts ',DOTS,' gives new transition state ',NTS+1
         NTS=NTS+1
         IF (NTS.GT.MAXTS) CALL TSDOUBLE
         ETS(NTS)=ETS(DOTS)
         FVIBTS(NTS)=FVIBTS(DOTS)
         IF (FRICTIONT) NEGEIG(NTS)=NEGEIG(DOTS)
         HORDERTS(NTS)=1          ! not valid in general ! see next line
         IF ((NATOMS.NE.7).OR.(.NOT.TWOD)) THEN
            PRINT*,'addperm assumes all minima have point group order 1 - quit'
            STOP
         ENDIF
         IXTS(NTS)=NIX
         IYTS(NTS)=NIY
         IZTS(NTS)=NIZ
         WRITE(UTS,REC=NTS) (LPOINTS(J2),J2=1,3*NATOMS)
         call flush(UTS,ISTAT)
C
C  Identify the connected minimum.
C
         DO J2=1,3*NATOMS
            LPOINTS(J2)=LPLUS(J2)
         ENDDO
         LDUMMY=LPOINTS(3*(NTAG-1)+1)
         LPOINTS(3*(NTAG-1)+1)=LPOINTS(3*(J1-1)+1)
         LPOINTS(3*(J1-1)+1)=LDUMMY
         LDUMMY=LPOINTS(3*(NTAG-1)+2)
         LPOINTS(3*(NTAG-1)+2)=LPOINTS(3*(J1-1)+2)
         LPOINTS(3*(J1-1)+2)=LDUMMY
         LDUMMY=LPOINTS(3*(NTAG-1)+3)
         LPOINTS(3*(NTAG-1)+3)=LPOINTS(3*(J1-1)+3)
         LPOINTS(3*(J1-1)+3)=LDUMMY

         CALL INERTIAWRAPPER(LPOINTS,NATOMS,angleAxis,NIX,NIY,NIZ)
         DO J2=1,NMIN
            IF ((ABS(EMIN(PLUS(DOTS))-EMIN(J2)).LT.EDIFFTOL).AND.
     1          (ABS(NIX-IXMIN(J2)).LT.IDIFFTOL).AND.
     2          (ABS(NIY-IYMIN(J2)).LT.IDIFFTOL).AND.
     3          (ABS(NIZ-IZMIN(J2)).LT.IDIFFTOL)) THEN
               WRITE(*,'(A,I5,A,I5)') 'permuting atoms 1 and ',J1,' for plus minimum gives old minimum ',J2
               PLUS(NTS)=J2
               IF (ENSEMBLE.EQ.'T') THEN
                  KPLUS(NTS)=LOG(1.0D0 * HORDERMIN(J2)  / (2.0D0 * PI*HORDERTS(NTS))) +
     1                       (FVIBMIN(J2)  - FVIBTS(NTS)) / 2.0D0 - (ETS(NTS) - EMIN(J2))/TEMPERATURE
                  IF (FRICTIONT) KPLUS(NTS)=KPLUS(NTS)+LOG(FRICTIONFAC(NEGEIG(NTS)))
               ELSE
                  IF (TEMPERATURE.GT.ETS(NTS)) THEN
                     KPLUS(NTS)  = LOG(1.0D0 * HORDERMIN(J2)  / (2*PI*HORDERTS(NTS))) +
     1               (FVIBMIN(J2)  - FVIBTS(NTS))/2 + (KAPPA-1)*LOG((TEMPERATURE-ETS(NTS))/(TEMPERATURE-EMIN(J2)))
                  ELSE
                      KPLUS(NTS)=-1.0D250
                  ENDIF
               ENDIF
               IF (ZSYM(1:2).EQ.'CA') KPLUS(NTS)=KPLUS(NTS)+30.66356D0
               IF (PLUS(NTS).EQ.MINUS(NTS)) KPLUS(NTS)=KPLUS(NTS)+LOG(2.0D0)  ! degenenerate rearrangement
!              IF (KSUM(J2).EQ.0.0D0) THEN
!                 IF (PLUS(NTS).NE.MINUS(NTS)) KSUM(J2)=KPLUS(NTS)
!              ELSE
!                 IF (PLUS(NTS).NE.MINUS(NTS)) KSUM(J2) =LOG(EXP(KSUM(J2)-KMEAN) + EXP(KPLUS(NTS) -KMEAN)) + KMEAN
!              ENDIF
               GOTO 11
            ENDIF
         ENDDO
         WRITE(*,'(A,I5,A,I5)') 'permuting atoms 1 and ',J1,' for plus minimum gives new minimum ',NMIN+1
         NMIN=NMIN+1
         IF (NMIN.GT.MAXMIN) CALL MINDOUBLE
         EMIN(NMIN)=EMIN(PLUS(DOTS))
         FVIBMIN(NMIN)=FVIBMIN(PLUS(DOTS))
         HORDERMIN(NMIN)=1          ! not valid in general !
         IXMIN(NMIN)=NIX
         IYMIN(NMIN)=NIY
         IZMIN(NMIN)=NIZ
         GPFOLD(NMIN)=0.0D0
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

         WRITE(UMIN,REC=NMIN) (LPOINTS(J2),J2=1,3*NATOMS)
         CALL FLUSH(UMIN,ISTAT)
         IF (CLOSEFILEST) OPEN(UNIT=UMINDATA,FILE='min.data',STATUS='UNKNOWN',POSITION='APPEND')
         WRITE(UMINDATA,'(2F20.10,I6,3F20.10)') EMIN(NMIN), FVIBMIN(NMIN), HORDERMIN(NMIN), IXMIN(NMIN), IYMIN(NMIN), IZMIN(NMIN)
         CALL FLUSH(UMINDATA,ISTAT)
         IF (CLOSEFILEST) CLOSE(UNIT=UMINDATA)

         PLUS(NTS)=NMIN
         IF (ENSEMBLE.EQ.'T') THEN
            KPLUS(NTS)=LOG(1.0D0 * HORDERMIN(NMIN)  / (2.0D0 * PI*HORDERTS(NTS))) +
     1                 (FVIBMIN(J2)  - FVIBTS(NTS)) / 2.0D0 - (ETS(NTS) - EMIN(NMIN))/TEMPERATURE
            IF (FRICTIONT) KPLUS(NTS)=KPLUS(NTS)+LOG(FRICTIONFAC(NEGEIG(NTS)))
         ELSE
            IF (TEMPERATURE.GT.ETS(NTS)) THEN
               KPLUS(NTS)  = LOG(1.0D0 * HORDERMIN(NMIN)  / (2*PI*HORDERTS(NTS))) +
     1            (FVIBMIN(NMIN)  - FVIBTS(NTS))/2 + (KAPPA-1)*LOG((TEMPERATURE-ETS(NTS))/(TEMPERATURE-EMIN(NMIN)))
            ELSE
               KPLUS(NTS)=-1.0D250
            ENDIF
         ENDIF
         IF (ZSYM(1:2).EQ.'CA') KPLUS(NTS)=KPLUS(NTS)+30.66356D0
         IF (PLUS(NTS).EQ.MINUS(NTS)) KPLUS(NTS)=KPLUS(NTS)+LOG(2.0D0)  ! degenenerate rearrangement
!        KSUM(NMIN)=0.0D0
!        IF (PLUS(NTS).NE.MINUS(NTS)) KSUM(NMIN)=KPLUS(NTS)
11       CONTINUE

         DO J2=1,3*NATOMS
            LPOINTS(J2)=LMINUS(J2)
         ENDDO
         LDUMMY=LPOINTS(3*(NTAG-1)+1)
         LPOINTS(3*(NTAG-1)+1)=LPOINTS(3*(J1-1)+1)
         LPOINTS(3*(J1-1)+1)=LDUMMY
         LDUMMY=LPOINTS(3*(NTAG-1)+2)
         LPOINTS(3*(NTAG-1)+2)=LPOINTS(3*(J1-1)+2)
         LPOINTS(3*(J1-1)+2)=LDUMMY
         LDUMMY=LPOINTS(3*(NTAG-1)+3)
         LPOINTS(3*(NTAG-1)+3)=LPOINTS(3*(J1-1)+3)
         LPOINTS(3*(J1-1)+3)=LDUMMY

         CALL INERTIAWRAPPER(LPOINTS,NATOMS,angleAxis,NIX,NIY,NIZ)
         DO J2=1,NMIN
            IF ((ABS(EMIN(MINUS(DOTS))-EMIN(J2)).LT.EDIFFTOL).AND.
     1          (ABS(NIX-IXMIN(J2)).LT.IDIFFTOL).AND.
     2          (ABS(NIY-IYMIN(J2)).LT.IDIFFTOL).AND.
     3          (ABS(NIZ-IZMIN(J2)).LT.IDIFFTOL)) THEN
               WRITE(*,'(A,I5,A,I5)') 'permuting atoms 1 and ',J1,' for minus minimum gives old minimum ',J2
               MINUS(NTS)=J2
               IF (ENSEMBLE.EQ.'T') THEN
                  KMINUS(NTS)=LOG(1.0D0 * HORDERMIN(J2) / (2.0D0 * PI*HORDERTS(NTS))) +
     1                     (FVIBMIN(J2) - FVIBTS(NTS)) / 2.0D0 - (ETS(NTS) - EMIN(J2))/TEMPERATURE
                  IF (FRICTIONT) KMINUS(NTS)=KMINUS(NTS)+LOG(FRICTIONFAC(NEGEIG(NTS)))
               ELSE
                  IF (TEMPERATURE.GT.ETS(NTS)) THEN
                     KMINUS(NTS) = LOG(1.0D0 * HORDERMIN(MINUS(NTS)) / (2*PI*HORDERTS(NTS))) +
     1                   (FVIBMIN(J2) - FVIBTS(NTS))/2 + 
     2                   (KAPPA-1)*LOG((TEMPERATURE-ETS(NTS))/(TEMPERATURE-EMIN(J2)))
                  ELSE
                     KMINUS(NTS)=-1.0D250
                  ENDIF
               ENDIF
               IF (ZSYM(1:2).EQ.'CA') KMINUS(NTS)=KMINUS(NTS)+30.66356D0
               IF (PLUS(NTS).EQ.MINUS(NTS)) KMINUS(NTS)=KMINUS(NTS)+LOG(2.0D0)  ! degenenerate rearrangement
!              IF (KSUM(J2).EQ.0.0D0) THEN
!                 IF (PLUS(NTS).NE.MINUS(NTS)) KSUM(J2)=KMINUS(NTS)
!              ELSE
!                 IF (PLUS(NTS).NE.MINUS(NTS)) KSUM(J2)=LOG(EXP(KSUM(J2)-KMEAN) + EXP(KMINUS(NTS)-KMEAN)) + KMEAN
!              ENDIF
               GOTO 12
            ENDIF
         ENDDO
         WRITE(*,'(A,I5,A,I5)') 'permuting atoms 1 and ',J1,' for minus minimum gives new minimum ',NMIN+1
         NMIN=NMIN+1
         IF (NMIN.GT.MAXMIN) CALL MINDOUBLE
         EMIN(NMIN)=EMIN(MINUS(DOTS))
         FVIBMIN(NMIN)=FVIBMIN(MINUS(DOTS))
         HORDERMIN(NMIN)=1          ! not valid in general !
         IXMIN(NMIN)=NIX
         IYMIN(NMIN)=NIY
         IZMIN(NMIN)=NIZ
         GPFOLD(NMIN)=0.0D0
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

         WRITE(UMIN,REC=NMIN) (LPOINTS(J2),J2=1,3*NATOMS)
         call flush(UMIN,ISTAT)
         IF (CLOSEFILEST) OPEN(UNIT=UMINDATA,FILE='min.data',STATUS='UNKNOWN',POSITION='APPEND')
         WRITE(UMINDATA,'(2F20.10,I6,3F20.10)') EMIN(NMIN), FVIBMIN(NMIN), HORDERMIN(NMIN), IXMIN(NMIN), IYMIN(NMIN),IZMIN(NMIN)
         CALL FLUSH(UMINDATA,ISTAT)
         IF (CLOSEFILEST) CLOSE(UNIT=UMINDATA)
         MINUS(NTS)=NMIN
         IF (ENSEMBLE.EQ.'T') THEN
            KMINUS(NTS)=LOG(1.0D0 * HORDERMIN(NMIN) / (2.0D0 * PI*HORDERTS(NTS))) +
     1               (FVIBMIN(NMIN) - FVIBTS(NTS)) / 2.0D0 - (ETS(NTS) - EMIN(NMIN))/TEMPERATURE
            IF (FRICTIONT) KMINUS(NTS)=KMINUS(NTS)+LOG(FRICTIONFAC(NEGEIG(NTS)))
         ELSE
            IF (TEMPERATURE.GT.ETS(NTS)) THEN
               KMINUS(NTS) = LOG(1.0D0 * HORDERMIN(NMIN) / (2*PI*HORDERTS(NTS))) +
     1             (FVIBMIN(NMIN) - FVIBTS(NTS))/2 + 
     2             (KAPPA-1)*LOG((TEMPERATURE-ETS(NTS))/(TEMPERATURE-EMIN(NMIN)))
            ELSE
               KMINUS(NTS)=-1.0D250
            ENDIF
         ENDIF
         IF (ZSYM(1:2).EQ.'CA') KMINUS(NTS)=KMINUS(NTS)+30.66356D0
         IF (PLUS(NTS).EQ.MINUS(NTS)) KMINUS(NTS)=KMINUS(NTS)+LOG(2.0D0)  ! degenenerate rearrangement
!        KSUM(NMIN)=0.0D0
!        IF (PLUS(NTS).NE.MINUS(NTS)) KSUM(NMIN)=KMINUS(NTS)
12       CONTINUE
C
C  Only now do we know the connected minima for sure.
C
         IF (CLOSEFILEST) OPEN(UNIT=UTSDATA,FILE='ts.data',STATUS='OLD',POSITION='APPEND')
         IF (IMFRQT) THEN
            WRITE(UTSDATA,'(2F20.10,3I10,4F20.10)') ETS(NTS),FVIBTS(NTS),HORDERTS(NTS),PLUS(NTS),MINUS(NTS),
     1                                          IXTS(NTS),IYTS(NTS),IZTS(NTS),NEGEIG(NTS)
         ELSE
            WRITE(UTSDATA,'(2F20.10,3I10,3F20.10)') ETS(NTS),FVIBTS(NTS),HORDERTS(NTS),PLUS(NTS),MINUS(NTS),
     1                                          IXTS(NTS),IYTS(NTS),IZTS(NTS)
         ENDIF
         CALL FLUSH(UTSDATA,ISTAT)
         IF (CLOSEFILEST) CLOSE(UNIT=UTSDATA)
C
C  Update ts pointers.
C
         POINTERP(NTS)=-1
         POINTERM(NTS)=-1
         TOPPOINTER(PLUS(NTS))=NTS
         TOPPOINTER(MINUS(NTS))=NTS

         DO J2=NTS-1,1,-1
            IF (PLUS(J2).EQ.PLUS(NTS)) THEN
               POINTERP(NTS)=J2
               GOTO 41
            ELSE IF (MINUS(J2).EQ.PLUS(NTS)) THEN
               POINTERP(NTS)=J2
               GOTO 41
            ENDIF
         ENDDO
41       CONTINUE

         DO J2=NTS-1,1,-1
            IF (PLUS(J2).EQ.MINUS(NTS)) THEN
               POINTERM(NTS)=J2
               GOTO 42
            ELSE IF (MINUS(J2).EQ.MINUS(NTS)) THEN
               POINTERM(NTS)=J2
               GOTO 42
            ENDIF
         ENDDO
42       CONTINUE
C
C  ts pointers have been updated.
C
10       CONTINUE
      ENDDO

      RETURN
      END
