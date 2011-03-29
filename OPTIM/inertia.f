C
C  GPL LICENSE INFO {{{
C 
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
C
C }}}
C 
C
      SUBROUTINE INERTIA(IT,Q)
! DOXYGEN {{{
!> \NAME INERTIA
!> \BRIEF BUILD INERTIA TENSOR IT FROM THE COORDINATES Q
!> \PARAM[OUT] IT
!> \PARAM[IN] Q
! }}}
      USE COMMONS
      USE KEY
      IMPLICIT NONE
      DOUBLE PRECISION IT(3,3),SCRATCH(9),Q(*), DIST
      INTEGER I, J, K, J1
C SUBROUTINE BODY {{{
C
C
      SCRATCH(1:9)=0.0D0
C      CMX=0.D0
C      CMY=0.D0
C      CMZ=0.D0
C      MOLWT=0.D0
C      DO I = 1,NATOMS
C         IF (TAGT.AND.(I.EQ.TAGNUM)) THEN
C            CMX = ATMASS(I)*Q(3*I-2)*TAGFAC+CMX
C            CMY = ATMASS(I)*Q(3*I-1)*TAGFAC+CMY
C            CMZ = ATMASS(I)*Q(3*I)*TAGFAC+CMZ
C            MOLWT = MOLWT+ATMASS(I)*TAGFAC
C         ELSE
C            CMX = ATMASS(I)*Q(3*I-2)+CMX
C            CMY = ATMASS(I)*Q(3*I-1)+CMY
C            CMZ = ATMASS(I)*Q(3*I)+CMZ
C            MOLWT = MOLWT+ATMASS(I)
C         ENDIF
CC        WRITE(*,'(A,I5,4F20.10)') 'I,ATMASS(I),=',I,ATMASS(I),Q(3*(I-1)+1),Q(3*(I-1)+2),Q(3*(I-1)+3)
C      ENDDO
C      SCRATCH(1) = CMX/MOLWT
C      SCRATCH(2) = CMY/MOLWT
C      SCRATCH(3) = CMZ/MOLWT
C      PRINT*,'MOLWT,SCRATCH=',MOLWT,SCRATCH(1),SCRATCH(2),SCRATCH(3)

      IT(1:3,1:3)=0.0D0
      IF (.NOT.MASST) THEN
         DO I=1,3
            KLOOP: DO K=1,NATOMS
               DO J1=1,NTAG
                  IF (K.EQ.TAGNUM(J1)) THEN
                     IT(I,I)=ATMASS(K)*(DIST(Q(3*K-2),SCRATCH(1))**2-Q(3*(K-1)+I)**2)*TAGFAC(J1)+IT(I,I)
                     CYCLE KLOOP
                  ENDIF
               ENDDO
               IT(I,I)=ATMASS(K)*(DIST(Q(3*K-2),SCRATCH(1))**2-Q(3*(K-1)+I)**2)+IT(I,I)
            ENDDO KLOOP
         ENDDO
         DO I=1,3
            DO J=1,3
               IF (I.NE.J) THEN
                  KLOOP2: DO K=1,NATOMS
                     DO J1=1,NTAG
                        IF (K.EQ.TAGNUM(J1)) THEN
                           IT(I,J)=-ATMASS(K)*(Q(3*(K-1)+I)*Q(3*(K-1)+J))*TAGFAC(J1)+IT(I,J)
                           CYCLE KLOOP2
                        ENDIF
                     ENDDO
                     IT(I,J)=-ATMASS(K)*(Q(3*(K-1)+I)*Q(3*(K-1)+J))+IT(I,J)
                  ENDDO KLOOP2
               ENDIF
            ENDDO
         ENDDO
      ELSE
         DO I=1,3
            KLOOP3: DO K=1,NATOMS
               DO J1=1,NTAG
                  IF (K.EQ.TAGNUM(J1)) THEN
                     IT(I,I)=(DIST(Q(3*K-2),SCRATCH(1))**2-Q(3*(K-1)+I)**2)*TAGFAC(J1)+IT(I,I)
                     CYCLE KLOOP3
                  ENDIF
               ENDDO
               IT(I,I)=(DIST(Q(3*K-2),SCRATCH(1))**2-Q(3*(K-1)+I)**2)+IT(I,I)
            ENDDO KLOOP3
         ENDDO
         DO I=1,3
            DO J=1,3
               IF (I.NE.J) THEN
                  KLOOP4: DO K=1,NATOMS
                     DO J1=1,NTAG
                        IF (K.EQ.TAGNUM(J1)) THEN
                           IT(I,J)=-(Q(3*(K-1)+I)*Q(3*(K-1)+J))*TAGFAC(J1)+IT(I,J)
                           CYCLE KLOOP4
                        ENDIF
                     ENDDO
                     IT(I,J)=-(Q(3*(K-1)+I)*Q(3*(K-1)+J))+IT(I,J)
                  ENDDO KLOOP4
               ENDIF
            ENDDO
         ENDDO
      ENDIF
      ! }}}
      RETURN
      END

      SUBROUTINE INERTIA2(DUMQ,ITX,ITY,ITZ)
      USE COMMONS
      USE KEY,ONLY : FREEZE, RBAAT
      IMPLICIT NONE
! SUBROUTINE PARAMETERS 
      DOUBLE PRECISION ITX,ITY,ITZ,DUMQ(3*NATOMS)
! LOCAL PARAMETERS 
      INTEGER J1, J2, J3, J4
      DOUBLE PRECISION IT(3,3), CMX, CMY, CMZ, VEC(3,3), MASST
C SUBROUTINE BODY {{{

      IF (RBAAT) NATOMS = NATOMS/2

      CMX=0.0D0
      CMY=0.0D0
      CMZ=0.0D0
      MASST=0.0D0
      DO J1=1,NATOMS
!        IF (TAGT.AND.(J1.EQ.TAGNUM)) THEN
            CMX=CMX+DUMQ(3*(J1-1)+1)*ATMASS(J1)*TAGFAC(J1)
            CMY=CMY+DUMQ(3*(J1-1)+2)*ATMASS(J1)*TAGFAC(J1)
            CMZ=CMZ+DUMQ(3*(J1-1)+3)*ATMASS(J1)*TAGFAC(J1)
            MASST=MASST+ATMASS(J1)*TAGFAC(J1)
!        ELSE
!           CMX=CMX+DUMQ(3*(J1-1)+1)*ATMASS(J1)
!           CMY=CMY+DUMQ(3*(J1-1)+2)*ATMASS(J1)
!           CMZ=CMZ+DUMQ(3*(J1-1)+3)*ATMASS(J1)
!           MASST=MASST+ATMASS(J1)
!        ENDIF
      ENDDO
      CMX=CMX/MASST
      CMY=CMY/MASST
      CMZ=CMZ/MASST
      IF (FREEZE) THEN
         CMX=0.0D0
         CMY=0.0D0
         CMZ=0.0D0
      ENDIF 
      DO J1=1,NATOMS
         DUMQ(3*(J1-1)+1)=DUMQ(3*(J1-1)+1)-CMX
         DUMQ(3*(J1-1)+2)=DUMQ(3*(J1-1)+2)-CMY
         DUMQ(3*(J1-1)+3)=DUMQ(3*(J1-1)+3)-CMZ
      ENDDO

      DO J1=1,3
         DO J2=1,3
            IT(J1,J2)=0.0D0
            J3LOOP1: DO J3=1,NATOMS
               DO J4=1,NTAG
                  IF (J3.EQ.TAGNUM(J4)) THEN
                     IT(J1,J2)=IT(J1,J2)-DUMQ(3*(J3-1)+J1)*DUMQ(3*(J3-1)+J2)*ATMASS(J3)*TAGFAC(J4)
                     CYCLE J3LOOP1
                  ENDIF
               ENDDO
               IT(J1,J2)=IT(J1,J2)-DUMQ(3*(J3-1)+J1)*DUMQ(3*(J3-1)+J2)*ATMASS(J3)
            ENDDO J3LOOP1
            IF (J1.EQ.J2) THEN
               J3LOOP2: DO J3=1,NATOMS
                  DO J4=1,NTAG
                     IF (J3.EQ.TAGNUM(J4)) THEN
                        IT(J1,J2)=IT(J1,J2)+(DUMQ(3*(J3-1)+1)**2+DUMQ(3*(J3-1)+2)**2+DUMQ(3*(J3-1)+3)**2)*ATMASS(J3)*TAGFAC(J4)
                        CYCLE J3LOOP2
                     ENDIF
                  ENDDO
                  IT(J1,J2)=IT(J1,J2)+(DUMQ(3*(J3-1)+1)**2+DUMQ(3*(J3-1)+2)**2+DUMQ(3*(J3-1)+3)**2)*ATMASS(J3)
               ENDDO J3LOOP2
            ENDIF
         ENDDO
      ENDDO

      CALL EIG(IT,VEC,3,3,0)
      ITX=IT(1,1)
      ITY=IT(2,2)
      ITZ=IT(3,3)

      IF (RBAAT) NATOMS = 2*NATOMS
C }}}
      RETURN
      END

