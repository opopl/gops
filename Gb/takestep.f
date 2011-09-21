
      SUBROUTINE TAKESTEP(NP)

      USE COMMONS

      IMPLICIT NONE

      DOUBLE PRECISION DPRAND, RANDOM, XMASS, YMASS, ZMASS, LOCALSTEP
      DOUBLE PRECISION :: DIST(3*NATOMS), DMAX, VMAX, VMIN, VMAX2, CMDIST(NATOMS), CMMAX

      DOUBLE PRECISION :: DUMMY
      DOUBLE PRECISION,PARAMETER :: PI=3.141592654D0 

      INTEGER J1, J2, JMAX, NP, J3, JMAX2
      
C  Calling CENTRE if NORESET is .TRUE. can lead to problems with COORDSO containing an atom
C  outside the permitted radius. Then it may be impossible to take a step that keeps all the
C  atoms inside.

! {{{
! COORDSO and VATO should now contain the correct COORDS and VAT even if
! these entries are changed by SYMMETRIZE or NEWRESET.
!
!     IF (.NOT.AMBERT) THEN
!        IF (ABS(COORDSO(1,NP)-COORDS(1,NP)).GT.1.0D-3) THEN
!           WRITE(MYUNIT,'(A,2G20.10)'),'takestep> WARNING - coordso will be changed: ',COORDSO(1,NP),COORDS(1,NP)
!        ENDIF
!        DO J1=1,3*(NATOMS-NSEED)
!           COORDSO(J1,NP)=COORDS(J1,NP)
!        ENDDO
!        DO J1=1,NATOMS
!           VATO(J1,NP)=VAT(J1,NP)
!        ENDDO
!     ENDIF
! }}}

C	csw34> FIND CENTRE OF COORDINATES (should be made a function!)

      XMASS=0.0D0; YMASS=0.0D0; ZMASS=0.0D0
      DO J1=1,NATOMS
         XMASS=XMASS+COORDS(3*(J1-1)+1,NP)
         YMASS=YMASS+COORDS(3*(J1-1)+2,NP)
         ZMASS=ZMASS+COORDS(3*(J1-1)+3,NP)
      ENDDO
      XMASS=XMASS/NATOMS; YMASS=YMASS/NATOMS; ZMASS=ZMASS/NATOMS

C  Find the most weakly bound atom, JMAX, the second most weakly bound atom, JMAX2,
C  and the pair energy of the most tightly bound atom, VMIN. An angular step is
C  taken for JMAX if its pair energy is > ASTEP*VMIN putting the atom at a radius of
C  DMAX (or CMMAX from CM of the cluster).

      DMAX=-1.0D0
      VMAX=-1.0D6
      VMAX2=-1.0D6
      VMIN=1.0D6
      CMMAX=-1.0D0
      DO J1=1,NATOMS
         J2=3*J1
         DIST(J1)= DSQRT( COORDS(J2-2,NP)**2+        COORDS(J2-1,NP)**2+        COORDS(J2,NP)**2)
         CMDIST(J1)=SQRT((COORDS(J2-2,NP)-XMASS)**2+(COORDS(J2-1,NP)-YMASS)**2+(COORDS(J2,NP)-ZMASS)**2)
         IF ((CMDIST(J1).GT.CMMAX).AND.(J1.LE.NATOMS-NCORE(NP))) CMMAX=CMDIST(J1)
         IF (DIST(J1).GT.DMAX) DMAX=DIST(J1)
            IF (VAT(J1,NP).GT.VMAX) THEN
               VMAX=VAT(J1,NP)
               JMAX=J1
            ELSE IF ((VAT(J1,NP).LT.VMAX).AND.(VAT(J1,NP).GT.VMAX2)) THEN
               VMAX2=VAT(J1,NP)
               JMAX2=J1
            ENDIF
         IF (VAT(J1,NP).LT.VMIN) VMIN=VAT(J1,NP)
      ENDDO
  

C	csw34> For the below loop, J1 is the atom number, and J2 is set to point to 3*J1, meaning
C	the z-coordinate of that atom. The y coordinate is then J2-1 and x is J2-2. 

      DO J1=1,NATOMS-NSEED
10       J2=3*J1
         LOCALSTEP=STEP(NP)

C  Tried changing to scale steps according to distance from CM. Maximum
C  allowed shift is linear with this distance. Worse.
C  Now try moving atoms according to how strongly bound they are.

           RANDOM=DPRAND()
           IF ((VMIN-VMAX.EQ.0.0D0).OR.(EFAC.EQ.0.0D0)) THEN
C	csw34> MAKE BASIC RANDOM CARTESIAN DISPLACEMENT MOVE
                 RANDOM=(DPRAND()-0.5D0)*2.0D0
                 COORDS(J2-2,NP)=COORDS(J2-2,NP)+LOCALSTEP*RANDOM
                 RANDOM=(DPRAND()-0.5D0)*2.0D0
                 COORDS(J2-1,NP)=COORDS(J2-1,NP)+LOCALSTEP*RANDOM
                 RANDOM=(DPRAND()-0.5D0)*2.0D0
                 COORDS(J2,NP)=COORDS(J2,NP)+LOCALSTEP*RANDOM
           ELSE 
              RANDOM=(DPRAND()-0.5D0)*2.0D0
              DUMMY=2.0D0+EAMP*TANH(-2.0D0*EFAC*(VAT(J1,NP)-(VMAX+VMIN)/2.0D0)/(VMAX-VMIN))
              COORDS(J2-2,NP)=COORDS(J2-2,NP)+LOCALSTEP*RANDOM*DUMMY
              RANDOM=(DPRAND()-0.5D0)*2.0D0
              COORDS(J2-1,NP)=COORDS(J2-1,NP)+LOCALSTEP*RANDOM*DUMMY
              RANDOM=(DPRAND()-0.5D0)*2.0D0
              COORDS(J2,NP)=COORDS(J2,NP)+LOCALSTEP*RANDOM*DUMMY
           ENDIF
               
      ENDDO

C  Preserve centre of mass if required.

      IF (CENT.AND.(.NOT.SEEDT)) CALL CENTRE2(COORDS(1:3*NATOMS,NP))
      IF (FIXCOM.AND.(.NOT.SEEDT)) CALL CENTRECOM(COORDS(:3*NATOMS,NP))
       
      RETURN
      END
