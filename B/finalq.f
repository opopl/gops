!op226>=================================== 
!op226> GPL License Info {{{ 
C   GMIN: A program for finding global minima
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of GMIN.
C
C   GMIN is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   GMIN is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
!op226>}}} 
!op226>=================================== 
      SUBROUTINE FINALQ
!op226>=================================== 
!op226> Declarations {{{ 
      USE COMMONS
      USE V
      IMPLICIT NONE

      INTEGER J1, J2, ITERATIONS, BRUN,QDONE, J3
      DOUBLE PRECISION POTEL, SCREENC(3*NATOMS), X(3*NATOMS), ENERGY, GRAD(3*NATOMS), TIME
      DOUBLE PRECISION SAVECSMNORM, DUMMY2, DIST2, RMAT(3,3), AVVAL, CSMGRAD(3), XTEMP(1:3*NATOMS)
      DOUBLE PRECISION DUMMY(3*NATOMS), AA(3)

      COMMON /MYPOT/ POTEL
!op226> End declarations }}} 
!op226>=================================== 
 
C
C  Make sure the lowest minima are tightly converged and then sort
C  them just to be on the safe side.
C
      CSMGUIDET=.FALSE.
      SHELLMOVES(1:NPAR)=.FALSE.
      IF (CUTT) CUTOFF=FINALCUTOFF
      SAVEQ=.FALSE.
      NQ(1)=0
      IF (FIELDT) FIELDT=.FALSE.
      IF (SEEDT) THEN
         SEEDT=.FALSE.
         NSEED=0
      ENDIF
      IF (SQUEEZET) SQUEEZET=.FALSE.
      MAXIT=MAXIT2
      DO J1=1,NSAVE
         IF (QMIN(J1).LT.1.0D10) THEN
            DO J2=1,3*NATOMS
               COORDS(J2,1)=QMINP(J1,J2)
            ENDDO
            NQ(1)=NQ(1)+1
            IF(UNFREEZEFINALQ) FROZEN(:)=.FALSE.  ! unfreeze everything before the final quenches
            CALL QUENCH(.TRUE.,1,ITERATIONS,TIME,BRUN,QDONE,SCREENC)
            WRITE(LFH,'(A,I6,A,F20.10,A,I5,A,F15.7,A,F12.2)') 'Final Quench ',NQ(1),' energy=',
     1                POTEL,' steps=',ITERATIONS,' RMS force=',RMS,' time=',TIME-TSTART

            QMIN(J1)=POTEL
            DO J2=1,3*NATOMS
               QMINP(J1,J2)=COORDS(J2,1)
            ENDDO
         ENDIF
      ENDDO

      NSAVE=NQ(1)
      CALL GSORT2(NSAVE,NATOMS)
C
C  Optionally sort the atoms from most bound to least bound according to VAT. {{{
C
       
      IF (SORTT) THEN
         DO J1=1,NSAVE
            IF (QMIN(J1).LT.0.0D0) THEN
               DO J2=1,3*NATOMS
                  COORDS(J2,1)=QMINP(J1,J2)
                  X(J2)=QMINP(J1,J2)
               ENDDO
               CALL POTENTIAL(X,GRAD,ENERGY,.FALSE.,.FALSE.)
               DO J2=1,NATOMS
                  VAT(J2,1)=VT(J2)
               ENDDO
               CALL SORT3(NATOMS,NATOMS,VAT(1:NATOMS,1),COORDS(1:3*NATOMS,1))
               DO J2=1,3*NATOMS
                  QMINP(J1,J2)=COORDS(J2,1)
               ENDDO
            ENDIF
         ENDDO
      ENDIF

!op226>}}} 

C     IF (DEBUG) THEN
         IF (TABOOT) THEN
            IF (NPAR.GT.1) THEN
               WRITE(LFH,'(A)') 'Taboo lists:'
               DO J1=1,NPAR
                  WRITE(LFH,'(A,G20.10)') 'Parallel run ',J1
                  WRITE(LFH,'(6F15.7)') (ESAVE(J2,J1),J2=1,NT(J1))
               ENDDO
            ELSE
               WRITE(LFH,'(A)') 'Taboo list:'
               WRITE(LFH,'(6F15.7)') (ESAVE(J2,1),J2=1,NT(1))
            ENDIF
         ENDIF
C     ENDIF
      RETURN
      END
