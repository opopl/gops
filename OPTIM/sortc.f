C   OPTIM: A PROGRAM FOR OPTIMIZING GEOMETRIES AND CALCULATING REACTION PATHWAYS
C   COPYRIGHT (C) 1999-2006 DAVID J. WALES
C   THIS FILE IS PART OF OPTIM.
C
C   OPTIM IS FREE SOFTWARE; YOU CAN REDISTRIBUTE IT AND/OR MODIFY
C   IT UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE AS PUBLISHED BY
C   THE FREE SOFTWARE FOUNDATION; EITHER VERSION 2 OF THE LICENSE, OR
C   (AT YOUR OPTION) ANY LATER VERSION.
C
C   OPTIM IS DISTRIBUTED IN THE HOPE THAT IT WILL BE USEFUL,
C   BUT WITHOUT ANY WARRANTY; WITHOUT EVEN THE IMPLIED WARRANTY OF
C   MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  SEE THE
C   GNU GENERAL PUBLIC LICENSE FOR MORE DETAILS.
C
C   YOU SHOULD HAVE RECEIVED A COPY OF THE GNU GENERAL PUBLIC LICENSE
C   ALONG WITH THIS PROGRAM; IF NOT, WRITE TO THE FREE SOFTWARE
C   FOUNDATION, INC., 59 TEMPLE PLACE, SUITE 330, BOSTON, MA  02111-1307  USA
C
      SUBROUTINE SORTC (NATMS, XX, ATNR, Y, NORD, X)
C     X IS A SCRATCH ARRAY OF THE SAME TYPE AND DIMENSION AS XX AND Y
C     ATNR ARE THE ATOMIX NUMBERS OF THE CENTERS (DUMMY ATOM = 0)
C
C     SORTS VECTOR OF NUCLEAR COORDINATES - TO CHECK FOR EQUIVALENCE
C     OF TWO ORIENTATIONS - NEEDS Q VECTOR AND ATOMIC NUMBER VECTOR (ATNR)
C
      IMPLICIT NONE
      INTEGER NATMS, ATNR(NATMS), I, J, JK, ILINE
      DOUBLE PRECISION X(3*NATMS),XX(3*NATMS),Y(3*NATMS),NORD(NATMS)
C
      ILINE(J)=1+J/3
C
C     SORT ON THE X - IF TWO X'S ARE EQUIVALENT, SORT ON Y AND SO ON.
C     
      DO I=1,3*NATMS
         X(I)=XX(I)
      ENDDO
C
C     FIRST GIVE DUMMY ATOMS RIDICULOUS SCRATCH COORDINATES - ENSURES
C     THAT THEY WILL WIND UP AT THE BOTTOM OF THE LIST
C
      DO 81 I=1,3*NATMS-2,3
         IF(ATNR(ILINE(I)) .EQ. 0)THEN
            DO J=0,2
               X(J+I) = -99995.
            ENDDO
         ENDIF
 81   CONTINUE
      JK=1
 429  J=1
      DO 96 I=1,3*NATMS-2,3
C
C     CONTINUE WITH SORTING.
C
         IF(X(I)-X(J).GT.1D-6)J=I
         IF(DABS(X(I)-X(J)).LT.1D-6)THEN
            IF(X(I+1)-X(J+1).GT.1D-6)J=I
            IF(DABS(X(I+1)-X(J+1)).LT.1D-6)THEN
               IF(X(I+2)-X(J+2).GT.1D-6)J=I
            ENDIF
         ENDIF
 96   CONTINUE
      DO    I=0,2
C
C     MASS-WEIGHT SORTED VECTOR - WILL ZERO ELEMENTS CORRESPONDING
C     TO DUMMY ATOMS SO THEY DONT MUCK UP THE SYMMETRY CHECKING.
C     
         Y(3*JK-2+I)=X(J+I)*ATNR(ILINE(J))
         X(J)=-99999.D0
      ENDDO
      NORD(JK)=(J+2)/3
      JK=JK+1
      IF(JK.EQ.NATMS+1)GO TO 999
      GO TO 429
 999  CONTINUE
      RETURN
      END
