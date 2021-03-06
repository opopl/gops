C   OPTIM: A program for optimizing geometries and calculating reaction pathways
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of OPTIM.
C
C   OPTIM is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   OPTIM is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
      SUBROUTINE UNRSETZSYMATMASS
      USE COMMONS
      USE MODUNRES
      IMPLICIT NONE

      INTEGER I1

      ALLOCATE(ATMASS(NATOMS))
      DO I1=1,NRES
         ZSYM(2*I1-1)='C'
         ZSYM(2*I1)='C'
         ATMASS(2*I1-1)=MASSES(1)
         ATMASS(2*I1)=MASSES(ITYPE(I1)+1)
      ENDDO
C JMC MAKE A MASSES FILE FOR PATHSAMPLE...
C     OPEN (UNIT=78,FILE='MASS',STATUS='UNKNOWN')
C     DO I1=1,NRES
C        WRITE(78,'(A3,F5.1)') ZSYM(2*I1-1),ATMASS(2*I1-1)
C        WRITE(78,'(A3,F5.1)') ZSYM(2*I1),ATMASS(2*I1)
C     END DO
C     CLOSE (UNIT=78)
 
      RETURN
      END


      SUBROUTINE UENERGY(X,GRAD,ETOT,GRADT,SECT)
      USE COMMONS
      USE MODUNRES
      IMPLICIT NONE

      DOUBLE PRECISION GRAD(3*NATOMS), X(3*NATOMS)
      REAL*8 EVDW,EVDW1,EVDW2,EES,EBE,ESC,ETORS,EHPB,EDIHCNSTR,ECORR,ETOT
      LOGICAL GRADT,SECT
      INTEGER ICALL
      COMMON /SRUTU/ ICALL

      NFL=0
      ICG=1
      ICALL=1

!CALL ZEROGRAD
!CALL ETOTAL(EVDW,EVDW1,EVDW2,EES,EBE,ESC,ETORS,EHPB,EDIHCNSTR,ECORR,ETOT)

      CALL INTGRAD(GRAD)
C JMC INTGRAD TAKES GRAD IN 'CARTS' (I.E. IN DC,DX: CA-CA AND SC-SC VECTORS) AND CONVERTS 
C IT TO A GRADIENT IN INTERNALS, WHICH IS PASSED BACK TO POTENTIAL ETC.

C JMC TESTING STUFF!!!
C     CALL ZEROGRAD
C     CALL CHECK_CARTGRAD
C     CALL CHECK_ECART
C     CALL CHECK_EINT

C     PRINT *,'UNRES ENERGY :'
C     PRINT *,ETOT

C     CALL CARTPRINT
C     CALL INT_FROM_CART(.TRUE.,.TRUE.)

      RETURN
      END SUBROUTINE UENERGY

      SUBROUTINE INTGRAD(GINT)
      USE COMMONS
      USE MODUNRES
      IMPLICIT NONE ! NOTE NOT PRESENT IN UNOPTIM.2.3

      DOUBLE PRECISION X(NINTS),G(NINTS)
      INTEGER UIPARM(1),ICALL,NF
      DOUBLE PRECISION URPARM(1)
      EXTERNAL FDUM
      DOUBLE PRECISION GINT(3*NATOMS)
      INTEGER J1
      COMMON /SRUTU/ ICALL

!CALL GEOM_TO_VAR(NVARU,X)
C     CALL VAR_TO_GEOM(NVARU,X)
C     CALL CHAINBUILD

      NF=1
      NFL=3
!CALL GRADIENT(NVARU,X,NF,G,UIPARM,URPARM,FDUM)

      GINT=0.0D0

      GINT(1:NVARU)=G(1:NVARU)

      RETURN
      END SUBROUTINE INTGRAD

      SUBROUTINE UPDATEDC
      USE COMMONS
      USE MODUNRES
      IMPLICIT NONE

      INTEGER J1

      DC=0.0D0
      DC_NORM=0.0D0

      DO J1=1,NRES-1
         DC(1,J1)=C(1,J1+1)-C(1,J1)
         DC(2,J1)=C(2,J1+1)-C(2,J1)
         DC(3,J1)=C(3,J1+1)-C(3,J1)
         DC_NORM(1,J1)=DC(1,J1)/VBL
         DC_NORM(2,J1)=DC(2,J1)/VBL
         DC_NORM(3,J1)=DC(3,J1)/VBL
C        PRINT *,MYDC(1,J1)/DC(1,J1),MYDC(2,J1)/DC(2,J1),MYDC(3,J1)/DC(3,J1)
C        PRINT *,MYDCNORM(1,J1)/DC_NORM(1,J1),MYDCNORM(2,J1)/DC_NORM(2,J1),MYDCNORM(3,J1)/DC_NORM(3,J1)
      END DO

      DO J1=NRES+1,2*NRES
         DC(1,J1)=C(1,J1)-C(1,J1-NRES)
         DC(2,J1)=C(2,J1)-C(2,J1-NRES)
         DC(3,J1)=C(3,J1)-C(3,J1-NRES)
         IF (ITYPE(J1-NRES).EQ.10) THEN   ! I.E. GLYCINE
            DC_NORM(1,J1)=DC(1,J1)
            DC_NORM(2,J1)=DC(2,J1)
            DC_NORM(3,J1)=DC(3,J1)
         ELSE
            DC_NORM(1,J1)=DC(1,J1)/DSC(ITYPE(J1-NRES))
            DC_NORM(2,J1)=DC(2,J1)/DSC(ITYPE(J1-NRES))
            DC_NORM(3,J1)=DC(3,J1)/DSC(ITYPE(J1-NRES))
         ENDIF
C        PRINT *,MYDC(1,J1)/DC(1,J1),MYDC(2,J1)/DC(2,J1),MYDC(3,J1)/DC(3,J1)
C        PRINT *,MYDCNORM(1,J1)/DC_NORM(1,J1),MYDCNORM(2,J1)/DC_NORM(2,J1),MYDCNORM(3,J1)/DC_NORM(3,J1)
      END DO
 
      RETURN
      END SUBROUTINE UPDATEDC
