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

! DJW STATIONARY POINT DUMP IN PATHSAMPLE FORMAT, ALTHOUGH WITH DUMMY ENTRIES
! FOR FREQUENCIES, POINT GROUP ORDERS AND MOMENTS OF INERTIA.

     SUBROUTINE DODUMPSP
     USE CONNECTDATA
     IMPLICIT NONE
     INTEGER NDUMMY, J1

! DUMP MINIMUM ENERGIES IN MIN.A.NEW, MIN.B.NEW AND MIN.DATA.NEW

     OPEN(24,FILE='MIN.A.NEW',STATUS='UNKNOWN')
     WRITE(24,'(2F20.10,I5,4F20.10)') MI(1)%DATA%E,0.0D0,1,1.0D0,1.0D0,1.0, 1.0D0
     CLOSE(24)
     OPEN(24,FILE='MIN.B.NEW',STATUS='UNKNOWN')
     WRITE(24,'(2F20.10,I5,4F20.10)') MI(2)%DATA%E,0.0D0,1,1.0D0,1.0D0,1.0, 1.0D0
     CLOSE(24)
     OPEN(24,FILE='MIN.DATA.NEW',STATUS='UNKNOWN')
     DO J1=3,NMIN
        WRITE(24,'(2F20.10,I5,4F20.10)') MI(J1)%DATA%E,0.0D0,1,1.0D0,1.0D0,1.0, 1.0D0
     ENDDO
     CLOSE(24)

! DUMP TS ENERGIES IN TS.DATA.NEW

     OPEN(24,FILE='TS.DATA.NEW',STATUS='UNKNOWN')
     DO J1=1,NTS
        WRITE(24,'(2F20.10,3I6,3F20.10)') TS(J1)%DATA%E,1.0D0,1,TS(J1)%DATA%P,TS(J1)%DATA%M,1.0D0,1.0D0,1.0D0
     ENDDO
     CLOSE(24)

! DUMP MINIMUM COORDINATES IN POINTS.MIN.NEW

     INQUIRE(IOLENGTH=NDUMMY) MI(1)%DATA%X(1:NOPT)
     OPEN(13,FILE='POINTS.MIN.NEW',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=NDUMMY)
     DO J1=1,NMIN
        WRITE(13,REC=J1) MI(J1)%DATA%X(1:NOPT)
     ENDDO
     CLOSE(13)

     OPEN(13,FILE='POINTS.TS.NEW',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=NDUMMY)
     DO J1=1,NTS
        WRITE(13,REC=J1) TS(J1)%DATA%X(1:NOPT)
     ENDDO
     CLOSE(13)
 
     END SUBROUTINE DODUMPSP

