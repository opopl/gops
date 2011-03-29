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
      SUBROUTINE KEGEN
      USE COMMONS, ONLY: GAUSSKK, GAUSSEE, GMODES, NATOMS, GKSMALL
      IMPLICIT NONE
      INTEGER IDUM, J, I
      DOUBLE PRECISION PI, GASDEV, DPRAND, SMALLEST

      IDUM = 1313417
      PI = ATAN(1.0D0)*4.0D0
      DO J=1,GMODES
!        EE(J) = 2*PI*RAN1(IDUM)  
         GAUSSEE(J) = 2.0D0*PI*DPRAND()  
         DO I=1,3*NATOMS
            GAUSSKK(I,J)=GASDEV(IDUM)
         ENDDO
      ENDDO
      DO I=1,3*NATOMS
         SMALLEST=1.0D100
         DO J=1,GMODES
            IF (DABS(GAUSSKK(I,J)).LT.SMALLEST) SMALLEST=GAUSSKK(I,J)
         ENDDO
         PRINT '(A,I8,G20.10)','I,SMALLEST KK=',I,SMALLEST
!
!   GKSMALL(I) SETS THE LENGTH SCALE FOR THE CORRESPONDING VARIABLE.
!   THE EFECTIVE RANGE IS 2*PI/GKSMALL(I)
!
         GKSMALL(I)=SMALLEST
      ENDDO

      RETURN
      END

      FUNCTION GASDEV(IDUM)
      IMPLICIT NONE
      INTEGER IDUM
      DOUBLE PRECISION GASDEV, DPRAND
      INTEGER ISET
      REAL FAC,GSET,RSQ,V1,V2,RAN1
      SAVE ISET,GSET
      DATA ISET/0/
      IF (ISET.EQ.0) THEN
! 1       V1=2.*RAN1(IDUM)-1.
!         V2=2.*RAN1(IDUM)-1.
1       V1=2.*DPRAND()-1.
        V2=2.*DPRAND()-1.
        RSQ=V1**2+V2**2
        IF(RSQ.GE.1..OR.RSQ.EQ.0.)GOTO 1
        FAC=SQRT(-2.*LOG(RSQ)/RSQ)
        GSET=V1*FAC
        GASDEV=V2*FAC
        ISET=1
      ELSE
        GASDEV=GSET
        ISET=0
      ENDIF
      RETURN
      END

      SUBROUTINE GFIELD(X,V,FUNCVALUE,GTEST)
      USE COMMONS, ONLY: NATOMS, GMODES, GAUSSKK, GAUSSEE
      IMPLICIT NONE
      LOGICAL GTEST,SECT
      INTEGER J, JJ, I
      DOUBLE PRECISION X(3*NATOMS), V(3*NATOMS), FUNCVALUE, NORM
      DOUBLE PRECISION SPROD(GMODES), SUMF

      NORM=SQRT(2.0D0)/SQRT(1.0D0*GMODES)

      DO J=1,GMODES
         SPROD(J)=0.0D0
         DO JJ=1,3*NATOMS
            SPROD(J)=SPROD(J)+GAUSSKK(JJ,J)*X(JJ)
         ENDDO
      ENDDO

      SUMF=0.0D0
      DO J=1,GMODES
          SUMF=SUMF+COS(SPROD(J)+GAUSSEE(J))
      ENDDO
 
      FUNCVALUE=NORM*SUMF

      DO JJ=1,3*NATOMS
         SUMF=0.0D0
         DO J=1,GMODES
            SUMF=SUMF-GAUSSKK(JJ,J)*SIN(SPROD(J)+GAUSSEE(J))
         ENDDO
         V(JJ)=SUMF*NORM
      ENDDO
  
      RETURN
       
      END
