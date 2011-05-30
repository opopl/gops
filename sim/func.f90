
MODULE GMIN.FUNC

IMPLICIT NONE
SAVE

INTEGER :: NUMBER_OF_ATOMS

CONTAINS

! Doxygen RAD  {{{
!
!> @name RAD
!> @brief Add energy and gradient correction terms for the container separately.
!> @param X
!> @param V
!> @param ENERGY
!> @param GTEST
!
! }}}
      SUBROUTINE RAD(X,V,ENERGY,GTEST)
      ! {{{

      USE COMMONS

      IMPLICIT NONE
      LOGICAL GTEST
      INTEGER J1, J3
      DOUBLE PRECISION X(3*NATOMS), DIST, V(3*NATOMS), 
     1                 ENERGY, DUMMYX, DUMMYY, DUMMYZ
      LOGICAL EVAP, EVAPREJECT
      COMMON /EV/ EVAP, EVAPREJECT

      IF (BSPT) RETURN ! container is accounted for in bspt by recounting previous configuration
      IF (PERIODIC) RETURN
      EVAP=.FALSE.
      EVAPREJECT=.FALSE.
      DO J1=1,NATOMS
         J3=3*J1
         DIST=X(J3-2)**2+X(J3-1)**2+X(J3)**2
         IF (DIST.GT.RADIUS) THEN
!           WRITE(MYUNIT,'(A,I5,5G20.10)') 'J1,DIST,RADIUS in rad = ',J1,DIST,RADIUS,X(J3-2),X(J3-1),X(J3)
            EVAP=.TRUE.
           ! IF (EVAP.AND.(BSWL.OR.BSPT)) then
           !     EVAPREJECT=.TRUE.
           !     IF (DEBUG) WRITE(MYUNIT,'(A,2G20.10)') 'EVAP: atom, radius=',J1,SQRT(DIST)
           !     RETURN
           ! ENDIF

            IF (DEBUG)  WRITE(MYUNIT,'(A,2G20.10,L10)') 'rad> EVAP: atom, radius, EVAP=',J1,SQRT(DIST),EVAP
C           PRINT*,'EVAP: atom, radius=',J1,SQRT(DIST)
CC           ENERGY=ENERGY+1.0D5*(DIST-RADIUS)**2
CC           IF (GTEST.AND.(.NOT.(SEEDT.AND.(J1.GT.NATOMS-NSEED).AND.FREEZECORE))) THEN
CC              DUMMYX=1.0D5*4.0D0*(DIST-RADIUS)*X(J3-2)
CC              DUMMYY=1.0D5*4.0D0*(DIST-RADIUS)*X(J3-1)
CC              DUMMYZ=1.0D5*4.0D0*(DIST-RADIUS)*X(J3)
CC              V(J3-2)=V(J3-2)+DUMMYX
CC              V(J3-1)=V(J3-1)+DUMMYY
CC              V(J3)=V(J3)+DUMMYZ
CC           ENDIF

             DIST=(SQRT(RADIUS)-0.5D0)/SQRT(DIST)
             X(J3-2)=X(J3-2)*DIST
             X(J3-1)=X(J3-1)*DIST
             X(J3)=X(J3)*DIST
!            WRITE(MYUNIT,'(A,3G20.10)') 'rad> reset coords: ',X(J3-2),X(J3-1),X(J3)
C
C  Put it back in at the opposite end of a diameter
C
C           X(J3-2)=-X(J3-2)*0.8D0
C           X(J3-1)=-X(J3-1)*0.8D0
C           X(J3)=-X(J3)*0.8D0
         ENDIF
      ENDDO
      ! }}}
      RETURN
      END

!  For rigid-body angle-axis coordinates, just move the fixed site

      SUBROUTINE RADR(X,V,ENERGY,GTEST)
      ! {{{
      USE COMMONS
      IMPLICIT NONE
      LOGICAL GTEST
      INTEGER J1, J3
      DOUBLE PRECISION X(3*NATOMS), DIST, V(3*NATOMS), 
     1                 ENERGY, DUMMYX, DUMMYY, DUMMYZ
      LOGICAL EVAP, evapreject
      COMMON /EV/ EVAP, evapreject

      IF (PERIODIC) RETURN
      EVAP=.FALSE.
      DO J1=1,NATOMS/2
         J3=3*J1
         DIST=X(J3-2)**2+X(J3-1)**2+X(J3)**2
C        WRITE(*,'(A,I6,5F15.5)') 'J1,DIST,coords,RADIUS in radr=',J1,DIST,X(J3-2),X(J3-1),X(J3),RADIUS
         IF (DIST.GT.RADIUS) THEN
            EVAP=.TRUE.
            WRITE(MYUNIT,'(A,I5,5G17.8)') 'EVAP: molecule, coords, dist, radius=',J1,X(J3-2),X(J3-1),X(J3),SQRT(DIST),SQRT(RADIUS)
C           IF (DEBUG) WRITE(*,'(A,I5,2G20.10)') 'EVAP: molecule, dist, radius=',J1,SQRT(DIST),SQRT(RADIUS)
            DIST=SQRT(RADIUS*0.9D0/DIST)
C           X(J3-2)=X(J3-2)*DIST
C           X(J3-1)=X(J3-1)*DIST
C           X(J3)=X(J3)*DIST
C
C  Put it back in at the opposite end of a diameter
C
            X(J3-2)=-X(J3-2)*0.8D0
            X(J3-1)=-X(J3-1)*0.8D0
            X(J3)=-X(J3)*0.8D0
         ENDIF
      ENDDO

      RETURN
      ! }}}
      END

SUBROUTINE OPENF(FILEHANDLE,MODE,FILENAME)
! {{{

INTEGER FILEHANDLE
CHARACTER (LEN=*) MODE,FILENAME

SELECTCASE(MODE)
        CASE(">")
                OPEN(FILEHANDLE,FILENAME,STATUS="UNKNOWN",FORM="FORMATTED")
        CASE("<")
        CASE(">>")
                OPEN(FILEHANDLE,FILENAME,STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
ENDSELECT
! }}}
ENDSUBROUTINE OPENF

SUBROUTINE MYRESET(NATOMS,NSEED)
! {{{
USE COMMONS,ONLY : COORDS,COORDSO,VAT,VATO

IMPLICIT NONE

INTEGER IA, NSEED, NATOMS

COORDS(1:NATOMS-NSEED)=COORDSO
VAT=VATO

RETURN
! }}}
END
 
SUBROUTINE TAKESTEP
      ! {{{

      USE COMMONS

      IMPLICIT NONE

      DO IA=1,NATOMS
         CALL GET_RND(RND,3,-1.0D0,1.0D0)
         COORDS(IA,1:3)=COORDS(IA,1:3)+STEP*RND(1:3)
      ENDDO
      
      RETURN
      ! }}}
END SUBROUTINE 

SUBROUTINE CENTRE2(R)
! {{{
USE COMMONS

IMPLICIT NONE

DOUBLE PRECISION RMASS(3), R(NATOMS,3)
INTEGER I,K

RMASS=SUM(R,DIM=1)/NATOMS

do K=1,3
R(:,K)=R(:,K)-RMASS(K)
ENDDO

IF (DEBUG) WRITE(LFH,'(A,3G20.10)') 'centre2> centre of mass reset to the origin from ',&                                                   (RMASS(I),I=1,3)
RETURN
! }}}
END


!> @name GSORT
!> @brief This subprogram performs a sort on the input data and
!> arranges it from smallest to biggest. The exchange-sort algorithm is used.
!
      SUBROUTINE GSORT(N,NATOMS)
! {{{

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NATOMS,N
      !DOUBLE PRECISION(:), INTENT(INOUT) ::
 
      INTEGER J1, L, J3, J2, NTEMP
      DOUBLE PRECISION TEMP, C
C
      DO 20 J1=1,N-1
         L=J1
         DO 10 J2=J1+1,N
            IF (QMIN(L).GT.QMIN(J2)) L=J2
10       CONTINUE
         TEMP=QMIN(L)
         QMIN(L)=QMIN(J1)
         QMIN(J1)=TEMP
         NTEMP=FF(L)
         FF(L)=FF(J1)
         FF(J1)=NTEMP
         DO J2=1,3*NATOMS
            C=QMINP(L,J2)
            QMINP(L,J2)=QMINP(J1,J2)
            QMINP(J1,J2)=C
         ENDDO
      ENDDO

      RETURN
! }}}
      END

SUBROUTINE GET_RND(RND,N,XMIN,XMAX)

IMPLICIT NONE

! random number vector
DOUBLE PRECISION :: RND(:),XMIN,XMAX,DX
! dimension of RND(:)
INTEGER N,I

DX=XMAX-XMIN

DO I=1,N 
        RND(I)=XMAX-DX*DPRAND()
ENDDO

RETURN

END 


ENDMODULE
