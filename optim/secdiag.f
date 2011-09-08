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
      SUBROUTINE SECDIAG(VEC,COORDS,ENERGY,GL,DIAG,GTEST,XRMS)
      USE COMMONS
      USE KEY
      USE modcharmm
      use porfuncs
      IMPLICIT NONE
      INTEGER J1
      LOGICAL GTEST, FPLUS, FMINUS
      DOUBLE PRECISION DUMMY3(3*NATOMS), DIFF, VEC(3*NATOMS), COORDS(3*NATOMS), DIAG2, DIAG3,
     1                 EPLUS,EMINUS,ENERGY,GRAD1(3*NATOMS),LOCALV(3*NATOMS),
     2                 RMS,DIAG,XRMS,GL(3*NATOMS),GRAD2(3*NATOMS), VECL, ZETA, PROJ

      DIFF=1.0D-3
      IF (AMHT) DIFF=1.0D-2
      IF (ZSYM(NATOMS).EQ.'GO') DIFF=2.0D-3
      IF (GAMESSUK.OR.GAMESSUS) DIFF=1.0D-2
      IF (CADPAC) DIFF=1.0D-2
      IF (CASTEP) DIFF=0.01D0
      IF (ONETEP) DIFF=0.01D0
      IF (CP2K) DIFF=0.001D0 
      IF (CPMD) DIFF=0.04D0
C     IF (CHRMMT) DIFF=5.0D-2
      IF (CHRMMT) DIFF=0.01D0
C     IF (DFTBT) DIFF=1.0D-3
C
C  Must read VEC into LOCALV because we are going to play with the vector in
C  question. This would mess up cases where we need to retry with a smaller
C  step size, because we cannot undo the previous step cleanly if the magnitude
C  of VEC is changed! 1/7/04 DJW
C
      LOCALV(1:NOPT)=VEC(1:NOPT)
C     PRINT '(A)','secdiag> vec before orthogopt'
C     PRINT '(6G20.10)',LOCALV(1:NOPT)
      IF (NFREEZE.LT.3) THEN
        CALL ORTHOGOPT(LOCALV,COORDS,.TRUE.)
      ELSE
         CALL VECNORM(LOCALV,NOPT)
      ENDIF
C     PRINT '(A)','secdiag> vec after orthogopt'
C     PRINT '(6G20.10)',LOCALV(1:NOPT)

      IF (FREEZE) THEN
         DO J1=1,NATOMS
            IF (FROZEN(J1)) THEN
               LOCALV(3*(J1-1)+1)=0.0D0
               LOCALV(3*(J1-1)+2)=0.0D0
               LOCALV(3*(J1-1)+3)=0.0D0
            ENDIF
         ENDDO
      ENDIF

      VECL=1.0D0
      ZETA=DIFF

      DO J1=1,NOPT
         DUMMY3(J1)=COORDS(J1)+ZETA*LOCALV(J1)
      ENDDO
      IF (CPMD) THEN
         INQUIRE(FILE='RESTART.1.plus',EXIST=FPLUS)
         IF (FPLUS) THEN
            CALL SYSTEM(' cp RESTART.1.plus RESTART.1 ')
         ELSE
            CALL SYSTEM(' cp RESTART.1.save RESTART.1 ')
         ENDIF
      ELSE IF (CASTEP) THEN
C        INQUIRE(FILE=SYS(1:LSYS) // '.wvfn.plus',EXIST=FPLUS)
C        IF (FPLUS) THEN
C           CALL SYSTEM(' cp ' // SYS(1:LSYS) // '.wvfn.plus ' // SYS(1:LSYS) // '.wvfn.1 ')
C        ELSE
C           CALL SYSTEM(' cp ' // SYS(1:LSYS) // '.wvfn.save ' // SYS(1:LSYS) // '.wvfn.1 ')
C        ENDIF
      ENDIF
C     PRINT*,'DUMMY3:'
      CALL POTENTIAL(DUMMY3,EPLUS,GRAD1,GTEST,.FALSE.,RMS,.FALSE.,.FALSE.)
      IF (CPMD) CALL SYSTEM(' cp RESTART.1 RESTART.1.plus ')
C     IF (CASTEP) CALL SYSTEM(' cp ' // SYS(1:LSYS) // '.wvfn.1 ' // SYS(1:LSYS) // '.wvfn.plus ')
      DO J1=1,NOPT
         DUMMY3(J1)=COORDS(J1)-ZETA*LOCALV(J1)
      ENDDO
      IF (CPMD) THEN
         INQUIRE(FILE='RESTART.1.minus',EXIST=FMINUS)
         IF (FMINUS) THEN
            CALL SYSTEM(' cp RESTART.1.minus RESTART.1 ')
         ELSE
            CALL SYSTEM(' cp RESTART.1.save RESTART.1 ')
         ENDIF
      ELSE IF (CASTEP) THEN
C        INQUIRE(FILE=SYS(1:LSYS) // '.wvfn.minus',EXIST=FMINUS)
C        IF (FMINUS) THEN
C           CALL SYSTEM(' cp ' // SYS(1:LSYS) // '.wvfn.minus ' // SYS(1:LSYS) // '.wvfn.1 ')
C        ELSE
C           CALL SYSTEM(' cp ' // SYS(1:LSYS) // '.wvfn.save ' // SYS(1:LSYS) // '.wvfn.1 ')
C        ENDIF
      ENDIF
C     PRINT*,'DUMMY3:'
      CALL POTENTIAL(DUMMY3,EMINUS,GRAD2,GTEST,.FALSE.,RMS,.FALSE.,.FALSE.)
      IF (CPMD) CALL SYSTEM(' cp RESTART.1 RESTART.1.minus ')
C     IF (CASTEP) CALL SYSTEM(' cp ' // SYS(1:LSYS) // '.wvfn.1 ' // SYS(1:LSYS) // '.wvfn.minus ')
      DIAG=(EPLUS+EMINUS-2.0D0*ENERGY)/(ZETA**2*VECL)
C     IF (DEBUG) WRITE(*,'(A,G15.5)') 'secdiag> DIFF',ZETA

      DIAG2=0.0D0
      DO J1=1,NOPT
         DIAG2=DIAG2+(GRAD1(J1)-GRAD2(J1))*LOCALV(J1)
C        WRITE(*,'(A,I4,4F20.10)') 'J1,GRAD1,GRAD2,LOCALV,DIAG2=',J1,GRAD1(J1),GRAD2(J1),LOCALV(J1),DIAG2
      ENDDO
      DIAG2=DIAG2/(2.0D0*ZETA)
      DIAG3=2*(DIAG-DIAG2/2)
C     IF (.NOT.GTEST) WRITE(*,'(A,6F20.10)') 'D,D2,D3,E+,E-,E=',DIAG,DIAG2,DIAG3,EPLUS,EMINUS,ENERGY
C
C  Although DIAG3 is a more accurate estimate of the diagonal second derivative, it
C  cannot be differentiated analytically.
C
      IF (GTEST) THEN
         DO J1=1,NOPT
            IF (NSECDIAG.EQ.2) THEN
               GL(J1)=(GRAD1(J1)-GRAD2(J1))/(ZETA*VECL**2)-2.0D0*DIAG2*LOCALV(J1)/VECL**2
            ELSE
               GL(J1)=(GRAD1(J1)-GRAD2(J1))/(ZETA*VECL**2)-2.0D0*DIAG*LOCALV(J1)/VECL**2
            ENDIF
!           WRITE(*,'(A,I4,4G16.7)') 'secdiag> J1,GRAD1,GRAD2,LOCALV,GL=',J1,GRAD1(J1),GRAD2(J1),LOCALV(J1),GL(J1)
         ENDDO
         IF (NFREEZE.LT.3) CALL ORTHOGOPT(GL,COORDS,.FALSE.)
C        PRINT *,'secdiag> before proj stuff GL:'
C        PRINT '(3F20.10)',GL(1:NOPT)
C        CALL ORTHOGOPT(GL,COORDS,.FALSE.) ! seems to do some good for MSEVB
C
C  Project out any component of the gradient along LOCALV (which is a unit vector)
C  This is a big improvement for DFTB.
C
         PROJ=0.0D0
         DO J1=1,NOPT
            PROJ=PROJ+GL(J1)*LOCALV(J1)
         ENDDO
         DO J1=1,NOPT
            GL(J1)=GL(J1)-PROJ*LOCALV(J1)
         ENDDO
         XRMS=0.0D0
         DO J1=1,NOPT
            XRMS=XRMS+GL(J1)**2
         ENDDO
         XRMS=DSQRT(XRMS/NOPT)
!        PRINT *,'secdiag> after proj stuff GL:'
!        PRINT '(3F20.10)',GL(1:NOPT)

         IF (DEBUG) WRITE(*,'(A,3G15.5,3G20.12,G10.3)') 'D,D2,D3,E+,E-,E,RMS=',DIAG,DIAG2,DIAG3,EPLUS,EMINUS,ENERGY,XRMS
         IF (DEBUG) WRITE(*,'(A,G20.10)') 'predicted gradient component=',(EPLUS-EMINUS)/(2*ZETA)
      ENDIF
!     PRINT '(A)','LOCALV:'
!     PRINT '(3G20.10)',LOCALV(1:NOPT)

      IF (NSECDIAG.EQ.2) DIAG=DIAG2

      RETURN
      END
