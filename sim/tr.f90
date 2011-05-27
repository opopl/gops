
!TAKESTEP

      ! the most weakly bound atom 
      INTEGER JMAX
      ! the second most weakly bound atom 
      INTEGER JMAX2
      ! pair energy of the most tightly bound atom 
      DOUBLE PRECISION VMIN
      ! for each particle, absolute distances
      DOUBLE PRECISION DIST(NATOMS)
      ! for each particle, distances relative to the centre of mass
      DOUBLE PRECISION CMDIST(NATOMS)

      FRPISQ = 4.D0*PI*PI
      NTRIESMAX=100

      RMASS=SUM(COORDS,dim=2)/NATOMS

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
      DO IA=1,NATOMS
         DIST(IA)= DSQRT(SUM(COORDS(IA,1:3)**2))
         CMDIST(IA)=SQRT(SUM(COORDS(IA,1:3)-RMASS(1:3))**2)
         IF ((CMDIST(IA).GT.CMMAX).AND.(IA.LE.NATOMS-NCORE)) CMMAX=CMDIST(IA)
         IF (DIST(IA).GT.DMAX) DMAX=DIST(IA)
            IF (VAT(IA).GT.VMAX) THEN
               VMAX=VAT(IA)
               JMAX=IA
            ELSE IF ((VAT(IA).LT.VMAX).AND.(VAT(IA).GT.VMAX2)) THEN
               VMAX2=VAT(IA)
               JMAX2=IA
            ENDIF
         IF (VAT(IA).LT.VMIN) VMIN=VAT(IA)
      ENDDO
! }}}


