
MODULE FUNC

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

      SUBROUTINE TAKESTEP
      ! {{{

      USE COMMONS

      IMPLICIT NONE

      ! the most weakly bound atom 
      INTEGER JMAX
      ! the second most weakly bound atom 
      INTEGER JMAX2
      ! pair energy of the most tightly bound atom 
      DOUBLE PRECISION VMIN
      ! for each particle, absolute distances
      DOUBLE PRECISION DIST(3*NATOMS)
      ! for each particle, distances relative to the centre of mass
      DOUBLE PRECISION CMDIST(NATOMS)

      FRPISQ = 4.D0*PI*PI
      NTRIESMAX=100

! Find centre of coordinates {{{

      XMASS=0.0D0; YMASS=0.0D0; ZMASS=0.0D0
      DO J1=1,NATOMS
         XMASS=XMASS+COORDS(3*(J1-1)+1,NP)
         YMASS=YMASS+COORDS(3*(J1-1)+2,NP)
         ZMASS=ZMASS+COORDS(3*(J1-1)+3,NP)
      ENDDO
      XMASS=XMASS/NATOMS; YMASS=YMASS/NATOMS; ZMASS=ZMASS/NATOMS

! }}}

!  Find JMAX, JMAX2, VMIN {{{
!
!  Find the most weakly bound atom, JMAX, the second most weakly bound atom, JMAX2,
!  and the pair energy of the most tightly bound atom, VMIN. An angular step is
!  taken for JMAX if its pair energy is > ASTEP*VMIN putting the atom at a radius of
!  DMAX (or CMMAX from CM of the cluster).
!

      DMAX=-1.0D0
      VMAX=-1.0D6
      VMAX2=-1.0D6
      VMIN=1.0D6
      CMMAX=-1.0D0
      DO J1=1,NATOMS
         J2=3*J1
         DIST(J1)= DSQRT( COORDS(J2-2)**2+COORDS(J2-1)**2+COORDS(J2)**2)
         CMDIST(J1)=SQRT((COORDS(J2-2)-XMASS)**2+(COORDS(J2-1)-YMASS)**2+(COORDS(J2)-ZMASS)**2)
         IF ((CMDIST(J1).GT.CMMAX).AND.(J1.LE.NATOMS-NCORE)) CMMAX=CMDIST(J1)
         IF (DIST(J1).GT.DMAX) DMAX=DIST(J1)
            IF (VAT(J1).GT.VMAX) THEN
               VMAX=VAT(J1)
               JMAX=J1
            ELSE IF ((VAT(J1).LT.VMAX).AND.(VAT(J1).GT.VMAX2)) THEN
               VMAX2=VAT(J1)
               JMAX2=J1
            ENDIF
         IF (VAT(J1).LT.VMIN) VMIN=VAT(J1)
      ENDDO
! }}}

      DO J1=1,NATOMS-NSEED
         IF ((COREFRAC.EQ.0.0D0).AND.(J1.GT.NATOMS-NCORE)) CYCLE ! no point taking a zero step
         NTRIES=0
10       J2=3*J1
         LOCALSTEP=STEP

              RANDOM=(DPRAND()-0.5D0)*2.0D0
              DUMMY=1.0D0+EAMP*TANH(-2.0D0*EFAC*(VAT(J1)-(VMAX+VMIN)/2.0D0)/(VMAX-VMIN))
              COORDS(J2-2)=COORDS(J2-2)+LOCALSTEP*RANDOM*DUMMY
              RANDOM=(DPRAND()-0.5D0)*2.0D0
              COORDS(J2-1)=COORDS(J2-1)+LOCALSTEP*RANDOM*DUMMY
              RANDOM=(DPRAND()-0.5D0)*2.0D0
              COORDS(J2)=COORDS(J2)+LOCALSTEP*RANDOM*DUMMY
      ENDDO

!  Preserve centre of mass if required.

      IF (CENT.AND.(.NOT.SEEDT)) CALL CENTRE2(COORDS(1:3*NATOMS))
      IF (FIXCOM.AND.(.NOT.SEEDT)) CALL CENTRECOM(COORDS(:3*NATOMS))
      
      RETURN
      ! }}}
      END SUBROUTINE 

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

ENDMODULE
