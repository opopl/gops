
!> @param[in]  X        dp(N)    input coordinates 
!> @param[out] GRADX     dp(N)    gradient 
!> @param[out] EREAL    dp         energy 
!> @param[out] RMS      dp         RMS 

      SUBROUTINE POTENTIAL(X,EREAL,GRADX,RMS,GRADT,SECT)
! declarations {{{

      IMPLICIT NONE

      ! subroutine parameters  {{{

      DOUBLE PRECISION, INTENT(OUT) :: EREAL, GRADX(:)
      DOUBLE PRECISION, INTENT(IN) :: X(:)
      LOGICAL GRADT, SECT
      DOUBLE PRECISION RMS
      ! }}}
      ! local parameters  {{{

      INTEGER NX,NR

      DOUBLE PRECISION,ALLOCATABLE :: R(:,:)
      DOUBLE PRECISION,ALLOCATABLE :: GRAD(:,:)
      ! }}}
      ! }}}

      NX=SIZE(X)
      NR=NX/3

      ALLOCATE(R(NR,3),GRAD(NR,3))

      R=RESHAPE(X,(/ NR,3 /))

      IF (BLNT .OR. PULLT) THEN 
            CALL EBLN(NR,R,ENERGY,GRAD,HESS,PTYPE,GRADT,SECT)
      ENDIF

      IF (PULLT) THEN
         EREAL=EREAL-PFORCE*R(PATOM1,3)+PFORCE*R(PATOM2,3)
         GRAD(PATOM1,3)=GRAD(PATOM1,3)-PFORCE
         GRAD(PATOM2,3)=GRAD(PATOM2,3)+PFORCE
      ENDIF

      X=RESHAPE(R,NX)
      GRADX=RESHAPE(GRAD,NX)
      RMS=MAX(DSQRT(SUM(GRADX**2)/NX),1.0D-100)

      RETURN
! }}}
      END
