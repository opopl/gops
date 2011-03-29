C   GMIN: A PROGRAM FOR FINDING GLOBAL MINIMA
C   COPYRIGHT (C) 1999-2006 DAVID J. WALES
C   THIS FILE IS PART OF GMIN.
C
C   GMIN IS FREE SOFTWARE; YOU CAN REDISTRIBUTE IT AND/OR MODIFY
C   IT UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE AS PUBLISHED BY
C   THE FREE SOFTWARE FOUNDATION; EITHER VERSION 2 OF THE LICENSE, OR
C   (AT YOUR OPTION) ANY LATER VERSION.
C
C   GMIN IS DISTRIBUTED IN THE HOPE THAT IT WILL BE USEFUL,
C   BUT WITHOUT ANY WARRANTY; WITHOUT EVEN THE IMPLIED WARRANTY OF
C   MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  SEE THE
C   GNU GENERAL PUBLIC LICENSE FOR MORE DETAILS.
C
C   YOU SHOULD HAVE RECEIVED A COPY OF THE GNU GENERAL PUBLIC LICENSE
C   ALONG WITH THIS PROGRAM; IF NOT, WRITE TO THE FREE SOFTWARE
C   FOUNDATION, INC., 59 TEMPLE PLACE, SUITE 330, BOSTON, MA  02111-1307  USA
C
C*************************************************************************
C
C  HERE WE CALCULATE THE EAMLJ POTENTIAL AND GRADIENT
C  (BASKES, PRL 27, 2592 (1999)). IT HAS THREE PARAMETERS A0, BETA AND Z0. 
C                                        
C*************************************************************************

      SUBROUTINE EAMLJ(XALL,V,ENERGY,GRADT)
      USE COMMONS
      IMPLICIT NONE 
      INTEGER NSIZE, I, J
      DOUBLE PRECISION XALL(3*NATOMS), ENERGY, X(NATOMS), Y(NATOMS), Z(NATOMS),
     +                 A0, BETA, Z0, RHOL(NATOMS), R(NATOMS,NATOMS), WW, V1, V2, 
     +                 V(3*NATOMS), FX(NATOMS), FY(NATOMS), FZ(NATOMS), BR, FUNC, 
     +                 RIJ, R6, WEXP, WLJ, WTOT 
      LOGICAL GRADT
      COMMON /EAMLJCOMM/ A0,BETA,Z0

      NSIZE=NATOMS

      DO I=1,NSIZE
        X(I)=XALL(3*(I-1)+1)
        Y(I)=XALL(3*(I-1)+2)
        Z(I)=XALL(3*(I-1)+3)
      ENDDO

C       CALCULATE ENERGY

      DO I=1,NSIZE
         DO J=I+1,NSIZE
            WW=(X(I)-X(J))**2+(Y(I)-Y(J))**2+(Z(I)-Z(J))**2
            WW=DSQRT(WW)
            R(I,J)=WW
            R(J,I)=R(I,J)
         ENDDO
      ENDDO
      
      V1=0.D0
      V2=0.D0

      DO I=1,NSIZE
         FX(I)=0.D0
         FY(I)=0.D0
         FZ(I)=0.D0
         RHOL(I)=0.D0
         DO J=1,NSIZE
            IF (I.NE.J) THEN
             RHOL(I)=RHOL(I)+DEXP(-BETA*(R(I,J)-1.D0))
            ENDIF
         ENDDO
         RHOL(I)=RHOL(I)/Z0
         V1=V1+FUNC(RHOL(I))

         DO J=I+1,NSIZE
            R6=1.D0/R(I,J)**6
            V2=V2+R6*(R6-2.D0)
            WW=DEXP(-BETA*(R(I,J)-1.D0))
            V2=V2-2.D0*FUNC(WW)/Z0
         ENDDO

      ENDDO

      ENERGY=V1+V2

C        PRINT*,ENERGY
 
        IF (.NOT.GRADT) RETURN

      DO I=1,NSIZE
         DO J=1,NSIZE
            IF (J.NE.I) THEN
             RIJ=R(J,I)

             BR=BETA*(RIJ-1.D0)
             WEXP=DEXP(-BR)

             WTOT=-0.5D0*(DLOG(RHOL(J)*RHOL(I)))
             WTOT=WTOT-BR
             WTOT=A0*BETA*WEXP*WTOT
             WLJ=12.D0/(RIJ**7)
             WLJ=WLJ*(1.D0-1.D0/RIJ**6)
             WTOT=WTOT+WLJ
C             WTOT=-WTOT
             FX(I)=FX(I)+(X(I)-X(J))*WTOT/RIJ
             FY(I)=FY(I)+(Y(I)-Y(J))*WTOT/RIJ
             FZ(I)=FZ(I)+(Z(I)-Z(J))*WTOT/RIJ
            ENDIF
         ENDDO
      ENDDO

        DO I=1,NSIZE
C          WRITE(*,'(3F20.10)') FX(I), FY(I), FZ(I)
          V(3*(I-1)+1)=FX(I)
          V(3*(I-1)+2)=FY(I)
          V(3*(I-1)+3)=FZ(I)
        ENDDO 

      RETURN
      END

C__________________________________________________________________________

      FUNCTION FUNC(X)
        IMPLICIT NONE
        DOUBLE PRECISION A0, BETA, Z0, X, FUNC
      COMMON/EAMLJCOMM/A0,BETA,Z0
      FUNC=A0*Z0*X
      FUNC=FUNC*(DLOG(X)-1.D0)
      FUNC=0.5D0*FUNC
      RETURN
      END
