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
MODULE TANGENT
     IMPLICIT NONE
     CONTAINS

     SUBROUTINE MEPTANGENT
          USE NEBDATA
          USE KEYNEB,ONLY: NIMAGE
          USE MODUNRES
          USE KEYTAU
          USE KEY, ONLY: BULKT, TWOD
          USE COMMONS, ONLY : PARAM1, PARAM2, PARAM3
          IMPLICIT NONE

          INTEGER :: J1,J2,K
          DOUBLE PRECISION :: WPLUS,WMINUS,DUMMY,TEMP1(NOPT),TEMP2(NOPT)
          DOUBLE PRECISION :: QINT(NINTS*(NIMAGE+2))
          DOUBLE PRECISION :: DIFFM(NINTS),DIFFP(NINTS)
          DOUBLE PRECISION :: PI=3.141592653589793D0

          IF ((TANTYPE.NE.1).AND.(BULKT)) THEN
             PRINT '(A)','meptangent> ERROR - minimum image distances only coded for TANTYPE 1'
             STOP
          ENDIF

          SELECT CASE (TANTYPE)
             CASE(1) ! Henkelman and Jonsson's improved tangent from JCP, 113, 9978, 2000
               DO J1=2,NIMAGE+1
                  CALL WS(EEE(J1-1),EEE(J1),EEE(J1+1),WMINUS,WPLUS)
                  TEMP1(1:NOPT)=0.0D0
                  TEMP2(1:NOPT)=0.0D0
                  IF (BULKT) THEN ! minimum image convention for distances!
                     DO K=1,NATOMS
                        TEMP1(3*(K-1)+1)=XYZ(NOPT*(J1-1)+3*(K-1)+1) - XYZ(NOPT*(J1-2)+3*(K-1)+1 ) &
   &                       -PARAM1*NINT((XYZ(NOPT*(J1-1)+3*(K-1)+1) - XYZ(NOPT*(J1-2)+3*(K-1)+1))/PARAM1)
                        TEMP1(3*(K-1)+2)=XYZ(NOPT*(J1-1)+3*(K-1)+2) - XYZ(NOPT*(J1-2)+3*(K-1)+2 ) &
   &                       -PARAM2*NINT((XYZ(NOPT*(J1-1)+3*(K-1)+2) - XYZ(NOPT*(J1-2)+3*(K-1)+2))/PARAM2)
                        IF (.NOT.TWOD) TEMP1(3*(K-1)+3)=XYZ(NOPT*(J1-1)+3*(K-1)+3) - XYZ(NOPT*(J1-2)+3*(K-1)+3 ) &
   &                       -PARAM3*NINT((XYZ(NOPT*(J1-1)+3*(K-1)+3) - XYZ(NOPT*(J1-2)+3*(K-1)+3))/PARAM3)

                        TEMP2(3*(K-1)+1)=XYZ(NOPT*J1+3*(K-1)+1) - XYZ(NOPT*(J1-1)+3*(K-1)+1 ) &
   &                       -PARAM1*NINT((XYZ(NOPT*J1+3*(K-1)+1) - XYZ(NOPT*(J1-1)+3*(K-1)+1))/PARAM1)
                        TEMP2(3*(K-1)+2)=XYZ(NOPT*J1+3*(K-1)+2) - XYZ(NOPT*(J1-1)+3*(K-1)+2 ) &
   &                       -PARAM2*NINT((XYZ(NOPT*J1+3*(K-1)+2) - XYZ(NOPT*(J1-1)+3*(K-1)+2))/PARAM2)
                        IF (.NOT.TWOD) TEMP2(3*(K-1)+3)=XYZ(NOPT*J1+3*(K-1)+3) - XYZ(NOPT*(J1-1)+3*(K-1)+3 ) &
   &                       -PARAM3*NINT((XYZ(NOPT*J1+3*(K-1)+3) - XYZ(NOPT*(J1-1)+3*(K-1)+3))/PARAM3)
                     ENDDO
                     TANVEC(1:NOPT,J1-1)=WMINUS*TEMP1(1:NOPT) + WPLUS*TEMP2(1:NOPT)
                  ELSE
                     TANVEC(:,J1-1) = WMINUS*( XYZ(NOPT*(J1-1)+1:NOPT*J1)         - XYZ(NOPT*(J1-2)+1:NOPT*(J1-1))         ) &
                                 &  + WPLUS *( XYZ(NOPT*J1+1:NOPT*(J1+1))         - XYZ(NOPT*(J1-1)+1:NOPT*J1    )         )
                  ENDIF
               ENDDO
             CASE(2) !NORMALIZED LINE SEGMENT (TANGENT ESTIMATED FROM TWO ADJACENT IMAGES)
               DO J1=1,NIMAGE
                    TANVEC(:,J1) = XYZ(NOPT*(J1+1)+1:NOPT*(J1+2)) - XYZ(NOPT*(J1-1)+1:NOPT*J1)
               ENDDO
             CASE(4) ! JMC FOR INTERNALS...
               TANVEC=0.0D0 ! INITIALISE TANVEC
               DO J1=1,NIMAGE+2
                  DO J2=1,NRES
                     C(1,J2)=XYZ(6*(J2-1)+1+NOPT*(J1-1))
                     C(2,J2)=XYZ(6*(J2-1)+2+NOPT*(J1-1))
                     C(3,J2)=XYZ(6*(J2-1)+3+NOPT*(J1-1))
                     C(1,J2+NRES)=XYZ(6*(J2-1)+4+NOPT*(J1-1))
                     C(2,J2+NRES)=XYZ(6*(J2-1)+5+NOPT*(J1-1))
                     C(3,J2+NRES)=XYZ(6*(J2-1)+6+NOPT*(J1-1))
                  ENDDO
                  CALL UPDATEDC
!CALL INT_FROM_CART(.TRUE.,.FALSE.)
!CALL GEOM_TO_VAR(NINTS,QINT(NINTS*(J1-1)+1:NINTS*J1))
               ENDDO
               DO J1=2,NIMAGE+1
                  CALL WS(EEE(J1-1),EEE(J1),EEE(J1+1),WMINUS,WPLUS)
                  DIFFM(:NINTS)=QINT(NINTS*(J1-1)+1:NINTS*J1) - QINT(NINTS*(J1-2)+1:NINTS*(J1-1))
                  DIFFP(:NINTS)=QINT(NINTS*J1+1:NINTS*(J1+1)) - QINT(NINTS*(J1-1)+1:NINTS*J1)
                  DO J2=1,NINTS
                     IF (DIFFM(J2).LT.-PI) DIFFM(J2) = DIFFM(J2)+2.0D0*PI
                     IF (DIFFP(J2).GT.PI) DIFFP(J2) = DIFFP(J2)-2.0D0*PI
                  ENDDO
                  TANVEC(:NINTS,J1-1) = WMINUS*DIFFM(:NINTS) + WPLUS*DIFFP(:NINTS)
               ENDDO
          END SELECT

          !  Normalise tangent vectors
          DO J1=1,NIMAGE
               IF (TANTYPE==4) THEN
                  DUMMY = SQRT( SUM(TANVEC(:NINTS,J1)**2) )
               ELSE
                  DUMMY = SQRT( SUM(TANVEC(:,J1)**2) )
               ENDIF

               IF (DUMMY==0.0D0) THEN
                    PRINT '(1x,a)', 'Tau norm is zero! Probably image density is too big...'
                    PRINT '(1x,a,i5)', ' This is image # ',j1
                    BADTAU=.TRUE.
                    RETURN
               ENDIF

               TANVEC(:,J1) = TANVEC(:,J1) / DUMMY
          ENDDO
     END SUBROUTINE MEPTANGENT

     SUBROUTINE WS(EL,EM,ER,WM,WP) ! ASSUMES THAT EL,EM AND ER ARE CHECKED FOR NAN ELSEWHERE ...
          IMPLICIT NONE
          DOUBLE PRECISION,INTENT(IN) :: EL,EM,ER
          DOUBLE PRECISION,INTENT(OUT) :: WM,WP

          ! El -->><<-- Em -->><<-- Er
          IF (EL < EM .AND. EM < ER) THEN
               WP=1.0D0;     WM=0.0D0
          ELSE IF (EL > EM .AND. EM > ER) THEN
               WP=0.0D0;     WM=1.0D0
          ELSE
               IF (EL < ER) THEN
                    WP=MAX(ABS(EL-EM),ABS(ER-EM));  WM=MIN(ABS(EL-EM),ABS(ER-EM))
               ELSE
                    WM=MAX(ABS(EL-EM),ABS(ER-EM));  WP=MIN(ABS(EL-EM),ABS(ER-EM))
               ENDIF
          ENDIF
     END SUBROUTINE WS

END MODULE TANGENT
