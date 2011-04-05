
MODULE FUNC

IMPLICIT NONE
SAVE

INTEGER :: NUMBER_OF_ATOMS

CONTAINS


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
