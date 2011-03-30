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
C
C*************************************************************************
C
C  Axilrod-Teller triple-dipole term.
C                                        
C*************************************************************************
C
      SUBROUTINE AXT(N,X,GRAD,POTEL,GRADT,ZSTAR)
      IMPLICIT NONE 
      INTEGER N, J1, J2, I, J, IJ
      DOUBLE PRECISION X(3*N), POTEL, P3, ZSTAR, DII, DJJ,
     1                 DIST(N,N), SDIST, DUMMY1, DUMMY2,
     2                 DOT(N,N), GRAD(3*N)
      LOGICAL GRADT
C
C  Distances.
C
      DO J1=1,N
         DIST(J1,J1)=0.0D0
         DO J2=J1+1,N
            SDIST=( X(3*(J2-1)+1)-X(3*(J1-1)+1) )**2 +
     1            ( X(3*(J2-1)+2)-X(3*(J1-1)+2) )**2 +
     2            ( X(3*(J2-1)+3)-X(3*(J1-1)+3) )**2 
            DIST(J1,J2)=1.0D0/DSQRT(SDIST)
            DIST(J2,J1)=DIST(J1,J2)
         ENDDO
      ENDDO
C
      DO J1=1,N
         DO J2=J1,N
            DOT(J2,J1)=X(3*(J2-1)+1)*X(3*(J1-1)+1)
     1                +X(3*(J2-1)+2)*X(3*(J1-1)+2)
     2                +X(3*(J2-1)+3)*X(3*(J1-1)+3)
            DOT(J1,J2)=DOT(J2,J1)
         ENDDO
      ENDDO
C
C  Calculate the energy
C
      P3=0.0D0
      DO I=1,N
         DII=DOT(I,I)
         DO J=I+1,N
            DJJ=DOT(J,J)
            DUMMY1=DOT(J,I)
            DUMMY2=DIST(I,J)
            DO IJ=J+1,N
               P3=P3+(1.0D0-3.0D0*
     1             (DUMMY1+DOT(IJ,J)-DOT(IJ,I)-DJJ)
     2            *(DUMMY1+DOT(IJ,IJ)-DOT(IJ,J)-DOT(IJ,I))
     3            *(DII+DOT(IJ,J)-DUMMY1-DOT(IJ,I))
     4            * (DUMMY2*DIST(I,IJ)*DIST(J,IJ))**2 ) 
     5            * (DUMMY2*DIST(I,IJ)*DIST(J,IJ))**3
            ENDDO
         ENDDO
      ENDDO
      P3=P3*ZSTAR
      POTEL=POTEL+P3
      IF (GRADT) CALL ATGRAD(N,X,GRAD,ZSTAR,DIST)

      RETURN
      END
C
C*************************************************************************
C
C  Here we calculate the analytic gradient for the Axilrod-Teller term.
C                                        
C*************************************************************************
C
      SUBROUTINE ATGRAD(N,X,V,ZSTAR,R2)
      IMPLICIT NONE 
      INTEGER N, J1, J2, J3, J4
      DOUBLE PRECISION X(3*N), ZSTAR, ABBC, ABAC, ACBC, ABI1, ABI2, ABI3, ACI1, ACI2, ACI3, BCI1, BCI2, BCI3,
     1                 R2(N,N), V(3*N), RX, RY, RZ, RR5(N,N),
     2                 VEC(3,N,N), DUMMY1, DUMMY2, DUMMY3, RR2(N,N), TDOT,
     3                 RAB, RRAB, RAC, RRAC, RBC, RRAB5, RRAC5, RRBC5
      DO J1=1,N
         R2(J1,J1)=0.0D0
         RR2(J1,J1)=0.0D0
         RR5(J1,J1)=0.0D0
         VEC(1,J1,J1)=0.0D0
         VEC(2,J1,J1)=0.0D0
         VEC(3,J1,J1)=0.0D0
         DO J2=J1+1,N
            R2(J2,J1)=1.0D0/R2(J2,J1)**2
            RR2(J2,J1)=1.0D0/R2(J2,J1)
            RR5(J2,J1)=1.0D0/R2(J2,J1)**(2.5D0)
            VEC(1,J2,J1)=X(3*(J2-1)+1)-X(3*(J1-1)+1)
            VEC(2,J2,J1)=X(3*(J2-1)+2)-X(3*(J1-1)+2)
            VEC(3,J2,J1)=X(3*(J2-1)+3)-X(3*(J1-1)+3)
            R2(J1,J2)=R2(J2,J1) 
            RR2(J1,J2)=RR2(J2,J1) 
            RR5(J1,J2)=RR5(J2,J1) 
            VEC(1,J1,J2)=-VEC(1,J2,J1)
            VEC(2,J1,J2)=-VEC(2,J2,J1)
            VEC(3,J1,J2)=-VEC(3,J2,J1)
         ENDDO
      ENDDO
C
C  Gradient
C
      DO J1=1,N
         DUMMY1=0.0D0
         DUMMY2=0.0D0
         DUMMY3=0.0D0
         DO J3=1,N
            RAB=R2(J3,J1)
            RRAB=RR2(J3,J1)
            RRAB5=RR5(J3,J1)
            ABI1=VEC(1,J3,J1)
            ABI2=VEC(2,J3,J1)
            ABI3=VEC(3,J3,J1)
            RX=VEC(1,J3,J1)
            RY=VEC(2,J3,J1)
            RZ=VEC(3,J3,J1)
            DO J4=J3+1,N
               ABAC=RX*VEC(1,J4,J1)+RY*VEC(2,J4,J1)+RZ*VEC(3,J4,J1)
               ABBC=RX*VEC(1,J4,J3)+RY*VEC(2,J4,J3)+RZ*VEC(3,J4,J3)
               ACBC=VEC(1,J4,J1)*VEC(1,J4,J3)+VEC(2,J4,J1)*VEC(2,J4,J3)+VEC(3,J4,J1)*VEC(3,J4,J3)
               TDOT=ABAC*ACBC*ABBC
               BCI1=VEC(1,J4,J3)
               BCI2=VEC(2,J4,J3)
               BCI3=VEC(3,J4,J3)
               ACI1=VEC(1,J4,J1)
               ACI2=VEC(2,J4,J1)
               ACI3=VEC(3,J4,J1)
               RBC=R2(J4,J3)
               RRBC5=RR5(J4,J3)
               RAC=R2(J4,J1)
               RRAC=RR2(J4,J1)
               RRAC5=RR5(J4,J1)*RRAB5*RRBC5
               DUMMY1=DUMMY1+ RRAC5 * ( 3*(ABAC*ACBC*BCI1 + ABBC*(ACBC*(ABI1 + ACI1) + ABAC*BCI1) + 
     1                  (ACI1*RAB + ABI1*RAC)*RBC) - 15*(ABI1*RRAB + ACI1*RRAC)*TDOT   )
               DUMMY2=DUMMY2+ RRAC5 * ( 3*(ABAC*ACBC*BCI2 + ABBC*(ACBC*(ABI2 + ACI2) + ABAC*BCI2) + 
     1                  (ACI2*RAB + ABI2*RAC)*RBC) - 15*(ABI2*RRAB + ACI2*RRAC)*TDOT   )
               DUMMY3=DUMMY3+ RRAC5 * ( 3*(ABAC*ACBC*BCI3 + ABBC*(ACBC*(ABI3 + ACI3) + ABAC*BCI3) + 
     1                  (ACI3*RAB + ABI3*RAC)*RBC) - 15*(ABI3*RRAB + ACI3*RRAC)*TDOT   )
            ENDDO
         ENDDO
         V(3*(J1-1)+1)=ZSTAR*DUMMY1+V(3*(J1-1)+1)
         V(3*(J1-1)+2)=ZSTAR*DUMMY2+V(3*(J1-1)+2)
         V(3*(J1-1)+3)=ZSTAR*DUMMY3+V(3*(J1-1)+3)
      ENDDO

      RETURN
      END
