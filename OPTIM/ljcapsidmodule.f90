MODULE LJCAPSIDMODULE

 INTEGER          :: I, J, J1, J2, J3, J4, J5, J6, K1, K2, OFFSET, REALNATOMS
 DOUBLE PRECISION, ALLOCATABLE :: RMIVEC(:,:,:), DPI1RMVEC(:,:,:), DPI2RMVEC(:,:,:), DPI3RMVEC(:,:,:)
 DOUBLE PRECISION, ALLOCATABLE :: PSCALEFAC1VEC(:),PSCALEFAC2VEC(:),EPSILON1(:,:,:),AEZR1(:,:,:), AEZR2(:,:,:)

 DOUBLE PRECISION :: ANGLE,ANGLE2,PI,TWOPI,SIGMA1,CUT,RON,RON2,RANGE2INV3

 DOUBLE PRECISION ::I3(3,3) 

END MODULE LJCAPSIDMODULE

SUBROUTINE INITIALISEPYGPERIODIC

USE KEY, ONLY: PARAMONOVPBCX,PARAMONOVPBCY,PARAMONOVPBCZ,PCUTOFF,PARAMONOVCUTOFF,&
                &       PYA1BIN,PYA2BIN,PYSIGNOT,PYEPSNOT,RADIFT,LJSITE,BLJSITE,PEPSILON1,&
                &       PSCALEFAC1,PSCALEFAC2,MAXINTERACTIONS,PYBINARYT,PYBINARYTYPE1,MYUNIT,VT,&
                &       PEPSILONATTR, PSIGMAATTR, LJSITEATTR, LJSITECOORDST, LJSITECOORDS, REALIGNXYZ
USE COMMONS, ONLY: NATOMS, RBSITE, NRBSITES
USE PYMODULE

IMPLICIT NONE
   WRITE(MYUNIT,*) ' INITIALISING VARIABLES FOR PY ',NATOMS,LJSITE,BLJSITE,PARAMONOVCUTOFF 
! ALLOCATE ARRAYS

    ALLOCATE(RMIVEC(NATOMS/2,3,3),DPI1RMVEC(NATOMS/2,3,3), DPI2RMVEC(NATOMS/2,3,3), DPI3RMVEC(NATOMS/2,3,3))
    ALLOCATE(PSCALEFAC1VEC(NATOMS/2),PSCALEFAC2VEC(NATOMS/2),EPSILON1(4,NATOMS/2,NATOMS/2))
    ALLOCATE(AEZR1(NATOMS/2,3,3), AEZR2(NATOMS/2,3,3))
    IF(.NOT.ALLOCATED(VT)) ALLOCATE(VT(NATOMS/2))
    IF(.NOT.ALLOCATED(PYA1BIN)) ALLOCATE(PYA1BIN(NATOMS/2,3))
    IF(.NOT.ALLOCATED(PYA2BIN)) ALLOCATE(PYA2BIN(NATOMS/2,3))


          VECSBF(1)=1.0D0
          VECSBF(2)=0.0D0
          VECSBF(3)=0.0D0

    IF(LJSITECOORDST) THEN
        VECSBF(:)=LJSITECOORDS(:)
        WRITE(MYUNIT,'(A,3F8.3)') ' REPULSIVE LJ SITE COORDINATES WILL BE ', LJSITECOORDS(:)
    END IF
         NRBSITES = 2
         ALLOCATE(RBSITE(NRBSITES,3))

         RBSITE(1,:)=0.0D0
         RBSITE(2,1)=LJSITECOORDS(1)
         RBSITE(2,2)=LJSITECOORDS(2)
         RBSITE(2,3)=LJSITECOORDS(3)

      I3(:,:)    = 0.D0
      AEZR1(:,:,:) = 0.D0
      AEZR2(:,:,:) = 0.D0

     PI=ATAN(1.0D0)*4.0
     TWOPI=2.0D0*PI

     IF(LJSITE) THEN
       EPSILON1(:,:,:)=PEPSILON1(1)   ! WE ARE GOING TO USE EPSILON1 FOR THE EXTRA LJ SITES
       SIGMA1(:)=1.0D0 ! PSIGMA1 IS NONEXISTENT FROM NOW ON, EXCEPT FOR THE ATTRACTIVE SECONDARY APEX SITES
       PSCALEFAC1VEC(:)=PSCALEFAC1(1)
       PSCALEFAC2VEC(:)=PSCALEFAC2(1)
       IF(LJSITEATTR) THEN !ATTRACTIVE SECONDARY APEX SITE IS TURNED ON
        ATTR(:)=1.0D0
        EPSILON1(1,:,:)=PEPSILONATTR(1)
        EPSILON1(2,:,:)=0.0D0!SQRT(PEPSILONATTR(1)*PEPSILONATTR(2))
        EPSILON1(3,:,:)=EPSILON1(2,:,:)
        EPSILON1(4,:,:)=PEPSILONATTR(2)

        SIGMA1(1)=PSIGMAATTR(1)
        SIGMA1(2)=0.5D0*(PSIGMAATTR(1)+PSIGMAATTR(2))
        SIGMA1(3)=SIGMA1(2)
        SIGMA1(4)=PSIGMAATTR(2)
       END IF
     ELSE
       EPSILON1(:,:,:)=0.0D0
       SIGMA1(:)=0.0D0
     END IF

! SANITY CHECKS
       IF(PYBINARYT.AND.LJSITE.AND..NOT.BLJSITE) THEN
        WRITE(MYUNIT,*) 'ERROR --- FOR BINARY PY SYSTEMS WITH EXTRA LJ SITES '// &
                        & 'YOU HAVE TO SPECIFY THE PARAMETERS FOR BOTH TYPES SEPARATELY (>3 ARGUMENTS AFTER EXTRALJSITE)! '
        STOP
       END IF
       IF(BLJSITE.AND..NOT.PYBINARYT) THEN
        WRITE(MYUNIT,*) 'ERROR --- BINARY LJ SITES SPECIFIED, BUT NO BINARY PY PARTICLES. '// &
                        & 'EXTRALJSITE SHOULD HAVE ONLY 3 ARGUMENTS!'
        STOP
       END IF


       IF(PYBINARYT.AND.BLJSITE) THEN
        DO J1=1,NATOMS/2
          IF(J1<=PYBINARYTYPE1) THEN
                PSCALEFAC1VEC(J1)=PSCALEFAC1(1)
                PSCALEFAC2VEC(J1)=PSCALEFAC2(1)
          ELSE
                PSCALEFAC1VEC(J1)=PSCALEFAC1(2)
                PSCALEFAC2VEC(J1)=PSCALEFAC2(2)
          END IF
        END DO

        DO J1=1,NATOMS/2-1
          DO J2=J1+1,NATOMS/2
           IF(J1<=PYBINARYTYPE1.AND.J2<=PYBINARYTYPE1) THEN
                EPSILON1(:,J1,J2)=PEPSILON1(1)
           ELSE IF(J1<=PYBINARYTYPE1.AND.J2>PYBINARYTYPE1) THEN
                EPSILON1(:,J1,J2)=PEPSILON1(3)
           ELSE
                EPSILON1(:,J1,J2)=PEPSILON1(2)
           END IF
          END DO
        END DO
!        WRITE(*,*) 'EPSILON1: ', EPSILON1(:,:,:)
       END IF
! CUTOFF STUFF FOR THE EXTRA LJ SITES
       CUT = PCUTOFF*PCUTOFF ! CUTOFF SQUARED
       RON = PCUTOFF*0.9D0
       RON2 = RON*RON
       RANGE2INV3=1.0D0/(CUT-RON2)**3
!     FROM INPUT PARAMETERS

      REALNATOMS = NATOMS/2
      OFFSET     = 3*REALNATOMS
      
      DO K1 = 1, 3
        I3(K1,K1) = 1.0D0
      ENDDO
      DO J1=1,REALNATOMS
       DO K1 = 1, 3
         AEZR1(J1,K1,K1) = 1.D0/(PYA1BIN(J1,K1)*PYA1BIN(J1,K1))
         AEZR2(J1,K1,K1) = 1.D0/(PYA2BIN(J1,K1)*PYA2BIN(J1,K1))
       ENDDO
      END DO
      IF(REALIGNXYZ) THEN
        CALL PYREALIGNXYZ
        STOP
      END IF
END SUBROUTINE INITIALISEPYGPERIODIC

! PY POTENTIAL, DC430'S IMPLEMENTATION
! WITH PBC AND CONTINUOUS CUTOFF ADDED
! PLUS EXTRA LJ SITE
      SUBROUTINE PYGPERIODIC (X, G, ENERGY, GTEST, STEST)

       USE KEY, ONLY: PARAMONOVPBCX,PARAMONOVPBCY,PARAMONOVPBCZ,PCUTOFF,PARAMONOVCUTOFF,&
                &       PYA1BIN,PYA2BIN,PYSIGNOT,PYEPSNOT,RADIFT,LJSITE,BLJSITE,PEPSILON1,&
                &       PSCALEFAC1,PSCALEFAC2,MAXINTERACTIONS,PYBINARYT,PYBINARYTYPE1,MYUNIT,VT,PEPSILONATTR,PSIGMAATTR
       USE COMMONS, ONLY: NATOMS
       USE PYMODULE
      IMPLICIT NONE

      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS)
      DOUBLE PRECISION :: RI(3), RJ(3), RIJ(3), NR(3), RIJSQ, ABSRIJ, P(3), THETA, THETA2, CT, ST
      DOUBLE PRECISION :: AE1(3,3), BE1(3,3), AE2(3,3), BE2(3,3), RM(3,3), APB(3,3), APBINV(3,3)
      DOUBLE PRECISION :: FCNT1, FCNT2, SRTFI1, SRTFI2, SRTFI(2), FMIN, LAMDAC1, LAMDAC2, ENERGY
      DOUBLE PRECISION :: RHO1, RHO1SQ, RHO16, RHO112, RHO2, RHO2SQ, RHO26
      DOUBLE PRECISION :: FCTR1, FCTR2, DVDF1, DVDF2 
      DOUBLE PRECISION :: DF1PI1, DF1PI2, DF1PI3, DF2PI1, DF2PI2, DF2PI3
      DOUBLE PRECISION :: DF1PJ1, DF1PJ2, DF1PJ3, DF2PJ1, DF2PJ2, DF2PJ3 
      DOUBLE PRECISION :: RMI(3,3), RMJ(3,3), E(3,3), ESQ(3,3), DE1(3,3), DE2(3,3), DE3(3,3)
      DOUBLE PRECISION :: DPI1RM(3,3), DPI2RM(3,3), DPI3RM(3,3), DPJ1RM(3,3), DPJ2RM(3,3), DPJ3RM(3,3)
      DOUBLE PRECISION :: DF1DR(3), DF2DR(3), DG1DR(3), DG2DR(3)
      DOUBLE PRECISION :: ARIBRJ(3), XC(3), XCMRI(3), XCMRJ(3), FIJ(3), TIJ(3), TJI(3)
      DOUBLE PRECISION :: D1ABEZ(3,3), D2ABEZ(3,3), D3ABEZ(3,3), D1ABE(3,3), D2ABE(3,3), D3ABE(3,3) 
      LOGICAL          :: GTEST


! SF344 ADDITIONS

      DOUBLE PRECISION :: CLJ1(2),CLJ2(2),CLJ3(2),CLJ4(2),CLJ5(2),CLJ6(2),CLJ7(2),CLJ8(2),CLJ11(2),CLJ12(2), &
                          & CLJ13(2),CLJ14(2),DFDR(2,3),DFP(2,6),DR(6),DCLJ1(2,12),DVDUMMY(12),LJ1(2),DUMMY,DUMMY1,DUMMY2,&
                          & DLJ1(2,12),VDUMMY
      DOUBLE PRECISION :: TERM2(MAXINTERACTIONS),TERM3(MAXINTERACTIONS), &
                          & XLJ(MAXINTERACTIONS,2,3),RLJVEC(MAXINTERACTIONS,3),RLJUNITVEC(MAXINTERACTIONS,3),&
                          & DRLJ(MAXINTERACTIONS,12),RLJ(MAXINTERACTIONS),RLJ2(MAXINTERACTIONS)
      DOUBLE PRECISION :: LLJ(12,MAXINTERACTIONS), DLLJ1(MAXINTERACTIONS,12)
      INTEGER          :: K
      LOGICAL          :: STEST
      DOUBLE PRECISION :: COSANGLE, SINANGLE, IMAT(3,3), TILDEMATRIX(3,3), ROTMAT(3,3), CROSSVECTOR(3),TEMPCRD(6)
     VT(1:NATOMS/2)=0.0D0
     TERM2(:)=1.0D0
     TERM3(:)=0.0D0
!       IF(PYBINARYT) THEN
!        DO J1=1,NATOMS/2-1
!          DO J2=J1+1,NATOMS/2
!           IF(J1<=PYBINARYTYPE1.AND.J2>PYBINARYTYPE1) THEN
!                EPSILON1(J1,J2)=PEPSILON1/10.0D0
!           ELSE
!                EPSILON1(J1,J2)=PEPSILON1
!          END DO
!        END DO
!       END IF

!      DO J1=1,REALNATOMS
!        J2=3*J1
!
!        IF (PARAMONOVPBCX) THEN
!!ENSURE X COMPONENT OF PARTICLE 1 VECTOR IS WITHIN BOXLX/2 OF ZERO. 
!!IF IT ISN'T THEN SUBTRACT INTEGER NUMBER OF BOXLX'S SUCH THAT IT IS.
!                X(J2-2)=X(J2-2)-BOXLX*NINT(X(J2-2)/BOXLX)
!        ENDIF
!
!        IF (PARAMONOVPBCY) THEN
!!ENSURE Y COMPONENT OF PARTICLE 1 VECTOR IS WITHIN BOXLY/2 OF ZERO. 
!!IF IT ISN'T THEN SUBTRACT INTEGER NUMBER OF BOXLY'S SUCH THAT IT IS.
!                X(J2-1)=X(J2-1)-BOXLY*NINT(X(J2-1)/BOXLY)
!        END IF
!
!        IF (PARAMONOVPBCZ) THEN
!!ENSURE Z COMPONENT OF PARTICLE 1 VECTOR IS WITHIN BOXLZ/2 OF ZERO. 
!!IF IT ISN'T THEN SUBTRACT INTEGER NUMBER OF BOXLZ'S SUCH THAT IT IS.
!                X(J2)=X(J2)-BOXLZ*NINT(X(J2)/BOXLZ)
!        ENDIF
!         
!      END DO

         ENERGY = 0.D0
         G(:)   = 0.D0
        
        DO J1=1, REALNATOMS

            J3      = 3*J1
            J5      = OFFSET + J3
            RI      = X(J3-2:J3)
            P       = X(J5-2:J5)
            ANGLE=SQRT(DOT_PRODUCT(P,P))
 
! TEMPORARY STUFF TO CREATE SPHEROIDAL STARTING GEOMETRY FOR GLOBAL OPTIMISATION
! STARTING COORDINATES ARE TAKEN FROM GLOBAL MINIMA OF THOMSON CLUSTERS OF THE APPROPRIATE SIZE
! - WE NEED TO FIGURE OUT THE ROTATION MATRIX THAT ROTATES THE RIGID BODY INTO THE VECTOR
! DEFINED BY THE POSITION OF THE PARTICLE (SINCE THAT IS A UNIT VECTOR)
! COS(ANGLE)=DOT_PRODUCT(VECSBF,COORDS(3*J-2:3*J))
!           CROSSVECTOR(:)=0.0D0
!           COSANGLE=DOT_PRODUCT(VECSBF,X(J3-2:J3))
!           CALL  CROSSOPT(VECSBF,X(J3-2:J3),CROSSVECTOR,0)
!           WRITE(*,*) 'CROSS VECTOR: ', CROSSVECTOR(:) 
!           SINANGLE=SQRT(DOT_PRODUCT(CROSSVECTOR,CROSSVECTOR))/(SQRT(DOT_PRODUCT(VECSBF,VECSBF))*SQRT(DOT_PRODUCT(X(J3-2:J3),X(J3-2:J3))))
!           IMAT(:,:)=0.0D0
!           WRITE(*,*) 'COSANGLE', COSANGLE, SINANGLE
!           WRITE(*,*) 'COSANGLE2+SINANGLE2=', COSANGLE*COSANGLE+SINANGLE*SINANGLE
!           DO I=1,3
!              IMAT(I,I)=1.0D0
!           END DO
!           TILDEMATRIX(1,1)=0.0D0
!           TILDEMATRIX(1,2)=-CROSSVECTOR(3)
!           TILDEMATRIX(1,3)= CROSSVECTOR(2)
!           TILDEMATRIX(2,1)=-TILDEMATRIX(1,2)
!           TILDEMATRIX(2,2)=0.0D0
!           TILDEMATRIX(2,3)=-CROSSVECTOR(1)
!           TILDEMATRIX(3,1)=-TILDEMATRIX(1,3)
!           TILDEMATRIX(3,2)=-TILDEMATRIX(3,2)
!           TILDEMATRIX(3,3)=0.0D0
!
!           ROTMAT=IMAT+SINANGLE*TILDEMATRIX+(1.0D0-COSANGLE)*MATMUL(TILDEMATRIX,TILDEMATRIX)
!
!           WRITE(*,'(A,3F8.3)') 'ORIGINAL VECTOR: ', VECSBF(:)
!           WRITE(*,'(A,3F8.3)') 'ORIGINAL VECTOR ROTATED ', MATMUL(VECSBF,ROTMAT)
!           WRITE(*,'(A,3F8.3)') 'SHOULD HAVE ROTATED INTO THIS: ', X(J3-2:J3)
!           TEMPCRD(1:6)=0.0D0
!           TEMPCRD(4)=1.0D0
!           CALL RBNEWROTGEOMMYORIENT(2,TEMPCRD,ROTMAT,0.0,0.0,0.0)
!           X(J5-2:J5)=TEMPCRD(4:6)
!           WRITE(*,'(A,3F8.3)') 'ANGLE-AXIS AFTER ROTATION: ',  X(J5-2:J5)
!           IF(J1==264) STOP



            IF(ANGLE>TWOPI) THEN
! NORMALISE ANGLE-AXIS COORDINATES
                X(J5-2:J5)=X(J5-2:J5)/ANGLE
                DO
                  ANGLE=ANGLE-TWOPI
                  IF(ANGLE<2*PI) EXIT
                END DO
! MULTIPLY WITH NEW ANGLE
                X(J5-2:J5)=X(J5-2:J5)*ANGLE
            END IF

            CALL RMDRVT(P, RMIVEC(J1,:,:), DPI1RMVEC(J1,:,:), DPI2RMVEC(J1,:,:), DPI3RMVEC(J1,:,:), .TRUE.)

        END DO

         DO J1 = 1, REALNATOMS - 1

            J3      = 3*J1
            J5      = OFFSET + J3
            RI      = X(J3-2:J3)
            P       = X(J5-2:J5)
!     ROTATION MATRIX

!            CALL RMDRVT(P, RMI, DPI1RM, DPI2RM, DPI3RM, GTEST)
            RMI(:,:)=RMIVEC(J1,:,:)
            DPI1RM(:,:)=DPI1RMVEC(J1,:,:)
            DPI2RM(:,:)=DPI2RMVEC(J1,:,:)
            DPI3RM(:,:)=DPI3RMVEC(J1,:,:)

            AE1 = MATMUL(RMI,(MATMUL(AEZR1(J1,:,:),(TRANSPOSE(RMI)))))

            IF (RADIFT) THEN

               AE2 = MATMUL(RMI,(MATMUL(AEZR2(J1,:,:),(TRANSPOSE(RMI)))))         

            ENDIF

!     BEGIN INNER LOOP OVER PARTICLES

            DO J2 = J1 + 1, REALNATOMS

               J4     = 3*J2
               J6     = OFFSET + J4
               RJ     = X(J4-2:J4) 
               P      = X(J6-2:J6)

!     ROTATION MATRIX

!               CALL RMDRVT(P, RMJ, DPJ1RM, DPJ2RM, DPJ3RM, GTEST)               
               RMJ(:,:)=RMIVEC(J2,:,:)
               DPJ1RM(:,:)=DPI1RMVEC(J2,:,:)
               DPJ2RM(:,:)=DPI2RMVEC(J2,:,:)
               DPJ3RM(:,:)=DPI3RMVEC(J2,:,:)
     
               BE1 = MATMUL(RMJ,(MATMUL(AEZR1(J2,:,:),(TRANSPOSE(RMJ)))))

               IF (RADIFT) THEN
   
                  BE2 = MATMUL(RMJ,(MATMUL(AEZR2(J2,:,:),(TRANSPOSE(RMJ)))))

               ENDIF

!     CALCULATE SEPARATION

               RIJ    = RI - RJ
               RIJSQ  = DOT_PRODUCT(RIJ,RIJ)
               ABSRIJ = DSQRT(RIJSQ)
               NR     = RIJ / ABSRIJ

                IF(PARAMONOVCUTOFF.AND.RIJSQ>CUT) GOTO 124

!     CALCULATE ECF

               CALL BRENTMIN (0.D0, 0.51D0, 1.D0, AE1, BE1, RIJ, LAMDAC1, FMIN)

               FCNT1   = - FMIN
               SRTFI1  = 1.D0 / DSQRT(FCNT1)
               APB     = LAMDAC1 * AE1 + (1.D0 - LAMDAC1) * BE1
               
               CALL MTRXIN (APB, APBINV)

               ARIBRJ =  LAMDAC1 * MATMUL(AE1,RI) + (1.D0 - LAMDAC1) * MATMUL(BE1,RJ)
               XC     =  MATMUL(APBINV, ARIBRJ)
               XCMRI  = XC - RI
               XCMRJ  = XC - RJ
               DF1DR  = - 2.D0 * LAMDAC1 * MATMUL(AE1,XCMRI)

               D1ABEZ = MATMUL(DPI1RM,AEZR1(J1,:,:))
               D1ABE  = MATMUL(D1ABEZ,(TRANSPOSE(RMI))) + MATMUL(RMI,(TRANSPOSE(D1ABEZ)))

               D2ABEZ = MATMUL(DPI2RM,AEZR1(J1,:,:))
               D2ABE  = MATMUL(D2ABEZ,(TRANSPOSE(RMI))) + MATMUL(RMI,(TRANSPOSE(D2ABEZ)))

               D3ABEZ = MATMUL(DPI3RM,AEZR1(J1,:,:))
               D3ABE  = MATMUL(D3ABEZ,(TRANSPOSE(RMI))) + MATMUL(RMI,(TRANSPOSE(D3ABEZ)))

               DF1PI1 = LAMDAC1*DOT_PRODUCT(XCMRI,MATMUL(D1ABE,XCMRI))
               DF1PI2 = LAMDAC1*DOT_PRODUCT(XCMRI,MATMUL(D2ABE,XCMRI))
               DF1PI3 = LAMDAC1*DOT_PRODUCT(XCMRI,MATMUL(D3ABE,XCMRI))

               D1ABEZ = MATMUL(DPJ1RM,AEZR1(J2,:,:))
               D1ABE  = MATMUL(D1ABEZ,(TRANSPOSE(RMJ))) + MATMUL(RMJ,(TRANSPOSE(D1ABEZ)))

               D2ABEZ = MATMUL(DPJ2RM,AEZR1(J2,:,:))
               D2ABE  = MATMUL(D2ABEZ,(TRANSPOSE(RMJ))) + MATMUL(RMJ,(TRANSPOSE(D2ABEZ)))

               D3ABEZ = MATMUL(DPJ3RM,AEZR1(J2,:,:))
               D3ABE  = MATMUL(D3ABEZ,(TRANSPOSE(RMJ))) + MATMUL(RMJ,(TRANSPOSE(D3ABEZ))) 
               
               DF1PJ1 = (1.D0-LAMDAC1)*DOT_PRODUCT(XCMRJ,MATMUL(D1ABE,XCMRJ))
               DF1PJ2 = (1.D0-LAMDAC1)*DOT_PRODUCT(XCMRJ,MATMUL(D2ABE,XCMRJ))
               DF1PJ3 = (1.D0-LAMDAC1)*DOT_PRODUCT(XCMRJ,MATMUL(D3ABE,XCMRJ))

               RHO1   = PYSIGNOT / (ABSRIJ - ABSRIJ*SRTFI1 + PYSIGNOT)
               RHO1SQ = RHO1*RHO1
               RHO16  = RHO1SQ*RHO1SQ*RHO1SQ
               RHO112 = RHO16 * RHO16

               FCTR1  = 0.5D0*ABSRIJ*SRTFI1/(FCNT1*PYSIGNOT)
               DG1DR  = (1.D0-SRTFI1)*NR/PYSIGNOT + FCTR1*DF1DR
               DVDF1  = -2.D0*RHO112*RHO1*FCTR1

               IF (RADIFT) THEN

                  CALL BRENTMIN (0.D0, 0.51D0, 1.D0, AE2, BE2, RIJ, LAMDAC2, FMIN)

                  FCNT2   = - FMIN
                  SRTFI2  = 1.D0 / DSQRT(FCNT2)
                  APB     = LAMDAC2 * AE2 + (1.D0 - LAMDAC2) * BE2

                  CALL MTRXIN (APB, APBINV)

                  ARIBRJ =  LAMDAC2 * MATMUL(AE2,RI) + (1.D0 - LAMDAC2) * MATMUL(BE2,RJ)
                  XC     =  MATMUL(APBINV, ARIBRJ)
                  XCMRI  = XC - RI
                  XCMRJ  = XC - RJ
                  DF2DR  = - 2.D0 * LAMDAC2 * MATMUL(AE2,XCMRI)

                  RHO2   = PYSIGNOT / (ABSRIJ - ABSRIJ*SRTFI2 + PYSIGNOT)
                  RHO2SQ = RHO2*RHO2
                  RHO26  = RHO2SQ*RHO2SQ*RHO2SQ
               
                  FCTR2  = 0.5D0*ABSRIJ*SRTFI2/(FCNT2*PYSIGNOT)
                  DG2DR  = (1.D0-SRTFI2)*NR/PYSIGNOT+FCTR2*DF2DR
                  DVDF2  = RHO26*RHO2*FCTR2

                  D1ABEZ = MATMUL(DPI1RM,AEZR2(J1,:,:))
                  D1ABE  = MATMUL(D1ABEZ,(TRANSPOSE(RMI))) + MATMUL(RMI,(TRANSPOSE(D1ABEZ)))

                  D2ABEZ = MATMUL(DPI2RM,AEZR2(J1,:,:))
                  D2ABE  = MATMUL(D2ABEZ,(TRANSPOSE(RMI))) + MATMUL(RMI,(TRANSPOSE(D2ABEZ)))

                  D3ABEZ = MATMUL(DPI3RM,AEZR2(J1,:,:))
                  D3ABE  = MATMUL(D3ABEZ,(TRANSPOSE(RMI))) + MATMUL(RMI,(TRANSPOSE(D3ABEZ)))

                  DF2PI1 = LAMDAC2*DOT_PRODUCT(MATMUL(XCMRI,D1ABE),XCMRI)
                  DF2PI2 = LAMDAC2*DOT_PRODUCT(MATMUL(XCMRI,D2ABE),XCMRI)
                  DF2PI3 = LAMDAC2*DOT_PRODUCT(MATMUL(XCMRI,D3ABE),XCMRI)

                  D1ABEZ = MATMUL(DPJ1RM,AEZR2(J2,:,:))
                  D1ABE  = MATMUL(D1ABEZ,(TRANSPOSE(RMJ))) + MATMUL(RMJ,(TRANSPOSE(D1ABEZ)))

                  D2ABEZ = MATMUL(DPJ2RM,AEZR2(J2,:,:))
                  D2ABE  = MATMUL(D2ABEZ,(TRANSPOSE(RMJ))) + MATMUL(RMJ,(TRANSPOSE(D2ABEZ)))

                  D3ABEZ = MATMUL(DPJ3RM,AEZR2(J2,:,:))
                  D3ABE  = MATMUL(D3ABEZ,(TRANSPOSE(RMJ))) + MATMUL(RMJ,(TRANSPOSE(D3ABEZ)))

                  DF2PJ1 = (1.D0-LAMDAC2)*DOT_PRODUCT(MATMUL(XCMRJ,D1ABE),XCMRJ)
                  DF2PJ2 = (1.D0-LAMDAC2)*DOT_PRODUCT(MATMUL(XCMRJ,D2ABE),XCMRJ)
                  DF2PJ3 = (1.D0-LAMDAC2)*DOT_PRODUCT(MATMUL(XCMRJ,D3ABE),XCMRJ)

               ELSE

                  RHO2   = RHO1
                  RHO26  = RHO16
                  DG2DR  = DG1DR
                  DVDF2  = RHO26*RHO2*FCTR1
                  DF2PI1 = DF1PI1
                  DF2PI2 = DF1PI2
                  DF2PI3 = DF1PI3
                  DF2PJ1 = DF1PJ1
                  DF2PJ2 = DF1PJ2
                  DF2PJ3 = DF1PJ3

                  SRTFI2 = SRTFI1
                  DF2DR  = DF1DR
                  RHO2SQ = RHO1SQ

               ENDIF


! CORRECTION TERMS TO THE POTENTIAL IF WE REQUIRE A CUTOFF AT RC. 
                 SRTFI(1)=SRTFI1
                 SRTFI(2)=SRTFI2
                 DFP(1,1)=DF1PI1
                 DFP(1,2)=DF1PI2
                 DFP(1,3)=DF1PI3
                 DFP(1,4)=DF1PJ1
                 DFP(1,5)=DF1PJ2
                 DFP(1,6)=DF1PJ3
                 DFP(2,1)=DF2PI1
                 DFP(2,2)=DF2PI2
                 DFP(2,3)=DF2PI3
                 DFP(2,4)=DF2PJ1
                 DFP(2,5)=DF2PJ2
                 DFP(2,6)=DF2PJ3
                 DFDR(1,:)=DF1DR(:)
                 DFDR(2,:)=DF2DR(:)
                 LJ1(1)=RHO1
                 LJ1(2)=RHO2
                   IF (PARAMONOVCUTOFF) THEN
! R/SQRT(F1(K)) = F(A,B), A PARAMETER ONLY DEPENDENT ON ORIENTATION!
                         DR(1)=1.0D0/ABSRIJ*(RI(1)-RJ(1))
                         DR(2)=1.0D0/ABSRIJ*(RI(2)-RJ(2))
                         DR(3)=1.0D0/ABSRIJ*(RI(3)-RJ(3))
                         DR(4)=-DR(1)
                         DR(5)=-DR(2)
                         DR(6)=-DR(3)
                     DO K=1,2
                          CLJ1(K)=PYSIGNOT/(PCUTOFF-ABSRIJ*SRTFI(K)+PYSIGNOT)
                          CLJ2(K)=CLJ1(K)**2
                          CLJ3(K)=CLJ2(K)*CLJ1(K)
                          CLJ4(K)=CLJ2(K)**2
                          CLJ5(K)=CLJ4(K)*CLJ1(K)
                          CLJ6(K)=CLJ4(K)*CLJ2(K)
                          CLJ7(K)=CLJ6(K)*CLJ1(K)
                          CLJ8(K)=CLJ6(K)*CLJ2(K)
                          CLJ11(K)=CLJ5(K)*CLJ6(K)
                          CLJ12(K)=CLJ6(K)**2
                          CLJ13(K)=CLJ12(K)*CLJ1(K)
                          CLJ14(K)=CLJ7(K)**2

                         DUMMY=CLJ1(K)/PYSIGNOT
                         DUMMY=DUMMY**2
                         DUMMY1=SRTFI(K)
                         DUMMY2=DUMMY1**3

                         DO J=1,3
                            DCLJ1(K,J) = 1.0D0*PYSIGNOT*DUMMY*(-1.0D0*DUMMY1*DR(J)+0.5D0*ABSRIJ*DUMMY2*DFDR(K,J))
                            DCLJ1(K,J+3)= -DCLJ1(K,J)
                         END DO
                         DO J=7,12
                            DCLJ1(K,J) =-1.0D0*PYSIGNOT*DUMMY*(0.5D0*ABSRIJ*DUMMY2*DFP(K,J-6)) !DERIVATIVES WRT TO ORIENTATION
                         END DO
                     END DO
                   END IF
 !     CALCULATE PY POTENTIAL ENERGY
123 CONTINUE 
        VDUMMY = 0.0D0
                 
               IF(LJSITE) THEN  
                DO K=1,MAXINTERACTIONS ! K=1 -- INTERACTION BETWEEN REPULSIVE PRIMARY 'APEX' SITES
                         ! K=2 AND K=3 -- INTERACTION BETWEEN SECONDARY AND PRIMARY 'APEX' SITES 
                         ! K=4 -- INTERACTION BETWEEN SECONDARY 'APEX' SITES (NORMAL LJ INTERACTION) 
                        ! TRYING TO MODIFY CODE TO ALLOW FOR BINARY SYSTEMS. 
                        ! APEX SITE HEIGHTS WILL BE DEFINED IN ABSOLUTE UNITS,
                        ! HENCE PYA1BIN(J1,1) ETC. WILL BE REMOVED FROM BELOW

                   IF(K==1) THEN
                        DUMMY1=PSCALEFAC1VEC(J1)!*PYA1BIN(J1,1)
                        DUMMY2=PSCALEFAC1VEC(J2)!*PYA1BIN(J2,1)
                   ELSE IF(K==2) THEN
                        DUMMY1=PSCALEFAC1VEC(J1)!*PYA1BIN(J1,1)
                        DUMMY2=-PSCALEFAC2VEC(J2)!*PYA1BIN(J2,1)
                   ELSE IF(K==3) THEN
                        DUMMY1=-PSCALEFAC2VEC(J1)!*PYA1BIN(J1,1)
                        DUMMY2=PSCALEFAC1VEC(J2)!*PYA1BIN(J2,1)
                   ELSE 
                        DUMMY1=-PSCALEFAC2VEC(J1)!*PYA1BIN(J1,1)
                        DUMMY2=-PSCALEFAC2VEC(J2)!*PYA1BIN(J2,1)
                   END IF
                        ! FIRST PARTICLE
                        XLJ(K,1,:)=RI+DUMMY1*MATMUL(RMI,VECSBF)    ! VECSBF: (1,0,0) IN THE BODY FRAME OF ELLIPSOID

                        ! SECOND PARTICLE
                        XLJ(K,2,:)=RJ+DUMMY2*MATMUL(RMJ,VECSBF)

                        ! SEPARATION BETWEEN THE LJ SITES
                        RLJ2(K)=(XLJ(K,2,1)-XLJ(K,1,1))**2+(XLJ(K,2,2)-XLJ(K,1,2))**2+(XLJ(K,2,3)-XLJ(K,1,3))**2
                        RLJ(K)=SQRT(RLJ2(K))
                        RLJVEC(K,1)=XLJ(K,2,1)-XLJ(K,1,1)
                        RLJVEC(K,2)=XLJ(K,2,2)-XLJ(K,1,2)
                        RLJVEC(K,3)=XLJ(K,2,3)-XLJ(K,1,3)

                        DUMMY=1.0D0/RLJ(K)
                        RLJUNITVEC(K,:)=RLJVEC(K,:)*DUMMY !/RLJ(K)

                        DRLJ(K,1)=DUMMY*(XLJ(K,2,1)-XLJ(K,1,1))         !DRLJ/DX1
                        DRLJ(K,2)=DUMMY*(XLJ(K,2,2)-XLJ(K,1,2))         !DRLJ/DY1
                        DRLJ(K,3)=DUMMY*(XLJ(K,2,3)-XLJ(K,1,3))         !DRLJ/DZ1
                        DRLJ(K,4)=-DRLJ(K,1)                               !DRLJ/DX2
                        DRLJ(K,5)=-DRLJ(K,2)                               !DRLJ/DY2
                        DRLJ(K,6)=-DRLJ(K,3)                               !DRLJ/DZ2
                        DRLJ(K,7) =-DUMMY*DUMMY1*DOT_PRODUCT(RLJVEC(K,:),MATMUL(DPI1RM,VECSBF)) !DRLJ/DPX1
                        DRLJ(K,8) =-DUMMY*DUMMY1*DOT_PRODUCT(RLJVEC(K,:),MATMUL(DPI2RM,VECSBF)) !DRLJ/DPY1
                        DRLJ(K,9) =-DUMMY*DUMMY1*DOT_PRODUCT(RLJVEC(K,:),MATMUL(DPI3RM,VECSBF)) !DRLJ/DPZ1
                        DRLJ(K,10) =  DUMMY*DUMMY2*DOT_PRODUCT(RLJVEC(K,:),MATMUL(DPJ1RM,VECSBF)) !DRLJ/DPX2
                        DRLJ(K,11) =  DUMMY*DUMMY2*DOT_PRODUCT(RLJVEC(K,:),MATMUL(DPJ2RM,VECSBF)) !DRLJ/DPY2
                        DRLJ(K,12) =  DUMMY*DUMMY2*DOT_PRODUCT(RLJVEC(K,:),MATMUL(DPJ3RM,VECSBF)) !DRLJ/DPZ2

              ! INTERACTION BETWEEN THE EXTRA LJ SITES:
                        LLJ(1,K)=SIGMA1(K)*DUMMY !/RLJ(K)
                        LLJ(2,K)=LLJ(1,K)**2
                        LLJ(3,K)=LLJ(2,K)*LLJ(1,K)
                        LLJ(4,K)=LLJ(2,K)**2
                        LLJ(5,K)=LLJ(4,K)*LLJ(1,K)
                        LLJ(6,K)=LLJ(4,K)*LLJ(2,K)
                        LLJ(7,K)=LLJ(6,K)*LLJ(1,K)
                        LLJ(11,K)=LLJ(5,K)*LLJ(6,K)
                        LLJ(12,K)=LLJ(6,K)*LLJ(6,K)

!                            DUMMY=1.0D0/RLJ(K)
!                            DUMMY=DUMMY**2
                            DO J=1,12
                                DLLJ1(K,J) =-SIGMA1(K)*DUMMY*DUMMY*DRLJ(K,J)
                            END DO

!                ADD CORRECTIONS TO
!                DERIVATIVES ARE ZERO AT RC, AND VANISH SMOOTHLY WITH NO DISCONTINUITIES. 
!                VDUMMY=EPSILON1*(LLJ12-0.0D0*LLJ6) !INNER HARD CORE
                 IF (PARAMONOVCUTOFF) THEN
                   IF (RLJ(K)>=PCUTOFF) THEN
!                           VDUMMY = 0.0D0
                           TERM2(K)=1.0D0
                           TERM3(K)=0.0D0
                   ELSE

                        ! WORK OUT THE SPLINE TERMS FOR SMOOTH CUTOFF OF EXTRA LJ SITES
                        IF(RLJ2(K)<RON2) THEN
                           TERM2(K)=1.0D0
                           TERM3(K)=0.0D0
                        ELSE IF(RLJ2(K)>RON2) THEN
                           TERM2(K)=(CUT-RLJ(K)**2)*(CUT-RLJ(K)**2)*(CUT+2.0D0*RLJ(K)**2-3.0D0*RON2)*RANGE2INV3
                           TERM3(K)=RLJ(K)*12.0D0*(CUT-RLJ(K)**2)*(RON2-RLJ(K)**2)*RANGE2INV3 ! D(TERM2)/DR
                        END IF
                   END IF
!                   VDUMMY = 0.0D0
                 END IF ! IF (PARAMONOVCUTOFF)
!                 VDUMMY=4.0D0*EPSILON1*(LLJ12-LLJ6)*TERM2 !EXTRA LJ SITE
                 VDUMMY=VDUMMY+4.0D0*EPSILON1(K,J1,J2)*TERM2(K)*(LLJ(12,K) - ATTR(K)*LLJ(6,K)) ! EXTRA LJ SITES
              END DO ! K=1,4
             END IF ! IF(LJSITE)

        IF(PARAMONOVCUTOFF) THEN
                !REPULSIVE POTENTIAL AND PERIODIC CUTOFF CORRECTIONS
                VDUMMY = VDUMMY + 4.0D0 * ( RHO112 + 6.0D0*CLJ12(1)*CLJ2(1)/RHO1SQ-7.0D0*CLJ12(1))
                    !                      LJ12(1)+(6.0D0*CLJ12(1))*CLJ2(1)/LJ2(1)-7.0D0*CLJ12(1)
                !ATTRACTIVE POTENTIAL AND PERIODIC CUTOFF CORRECTIONS
                VDUMMY = VDUMMY + 4.0D0 * (- RHO26 - 3.0D0* CLJ6(2)*CLJ2(2)/RHO2SQ+4.0D0* CLJ6(2))
                    !                      (-LJ6(2)+(-3.0D0*CLJ6(2))*CLJ2(2)/LJ2(2) +4.0D0*CLJ6(2))
        ELSE
                VDUMMY = VDUMMY + 4.0D0 * (RHO112 - RHO26)
        END IF

               ENERGY = ENERGY + VDUMMY
               VT(J1) = VT(J1) + VDUMMY
               VT(J2) = VT(J2) + VDUMMY        ! PAIR POTENTIALS

            DVDUMMY(:) = 0.0D0

!     CALCULATE GRADIENT
             FIJ = 0.0D0
             TIJ = 0.0D0
             TJI = 0.0D0

          IF(PARAMONOVCUTOFF) THEN
            DVDUMMY(:) = 0.0D0
            ! WITH RESPECT TO CARTESIANS
            DO J=1,3
              IF(LJSITE) THEN
               DO K=1,MAXINTERACTIONS
!               DVDUMMY(J) = 24.0D0*EPSILON1*(2.0D0*LLJ11*DLLJ1(J)-LLJ5*DLLJ1(J))*1.0D0*TERM2 + &
!                          & 4.0D0*EPSILON1*(LLJ12-LLJ6)*TERM3*DRLJ(J)! EXTRA LJ SITE DERIVATIVES
                DVDUMMY(J) = DVDUMMY(J) + 4.0D0*EPSILON1(K,J1,J2)*(12.0D0*LLJ(11,K)*DLLJ1(K,J))*TERM2(K) + &
                          & 4.0D0*EPSILON1(K,J1,J2)*LLJ(12,K)*TERM3(K)*DRLJ(K,J)! EXTRA LJ SITE DERIVATIVES, NOW ONLY REPULSIVE
                DVDUMMY(J) = DVDUMMY(J) - ATTR(K)*(4.0D0*EPSILON1(K,J1,J2)*(6.0D0*LLJ(5,K)*DLLJ1(K,J))*TERM2(K) + &
                          & 4.0D0*EPSILON1(K,J1,J2)*LLJ(6,K)*TERM3(K)*DRLJ(K,J)) ! ATTRACTIVE SECONDARY APEX SITE
               END DO
              ELSE
               DVDUMMY(J) = 0.0D0
              END IF
!                DVDUMMY(J) = DLLJ1(J)
               DVDUMMY(J) = DVDUMMY(J) + 24.0D0 * ( 2.0D0*RHO1SQ*DG1DR(J)*(RHO16*RHO1SQ*RHO1SQ*RHO1-CLJ14(1)/(RHO1SQ*RHO1))&
                   &+14.0D0*DCLJ1(1,J)*(CLJ13(1)/RHO1SQ-CLJ11(1))) !REPULSIVE DERIVATIVE

               DVDUMMY(J) = DVDUMMY(J) + 24.0D0 * (-1.0D0*RHO2SQ*DG2DR(J)*(RHO2SQ*RHO2SQ*RHO2-CLJ8(2)/(RHO2SQ*RHO2))&
                   &-4.0D0*DCLJ1(2,J)*(CLJ7(2)/RHO2SQ-CLJ5(2))) !ATTRACTIVE DERIVATIVE
!               DVDUMMY(J) = DCLJ1(1,J)
               FIJ(J) = DVDUMMY(J)
            END DO

            DO K=1,2
                DUMMY=LJ1(K)/PYSIGNOT
                DUMMY=DUMMY**2
                DUMMY1=SRTFI(K)
                DUMMY2=DUMMY1**3
               DO J=7,12
                  DLJ1(K,J) =-1.0D0*PYSIGNOT*DUMMY*(0.5D0*ABSRIJ*DUMMY2*DFP(K,J-6))
               END DO
            END DO
            ! WITH RESPECT TO ORIENTATION VECTORS
            DO J=7,12
             IF(LJSITE) THEN
              DO K=1,MAXINTERACTIONS
!               DVDUMMY(J) = 24.0D0*EPSILON1*(2.0D0*LLJ1**11*DLLJ1(J)-LLJ1**5*DLLJ1(J))*TERM2 + &
!                               & 4.0D0*EPSILON1*(LLJ12-LLJ6)*TERM3*DRLJ(J)! EXTRA LJ SITE DERIVATIVES
               DVDUMMY(J) = DVDUMMY(J) + 4.0D0*EPSILON1(K,J1,J2)*(12.0D0*LLJ(11,K)*DLLJ1(K,J))*TERM2(K) + &
                               & 4.0D0*EPSILON1(K,J1,J2)*LLJ(12,K)*TERM3(K)*DRLJ(K,J)! EXTRA LJ SITE DERIVATIVES, NOW ONLY REPULSIVE
               DVDUMMY(J) = DVDUMMY(J) - ATTR(K)*(4.0D0*EPSILON1(K,J1,J2)*(6.0D0*LLJ(5,K)*DLLJ1(K,J))*TERM2(K) + &
                          & 4.0D0*EPSILON1(K,J1,J2)*LLJ(6,K)*TERM3(K)*DRLJ(K,J)) ! ATTRACTIVE SECONDARY APEX SITE
              END DO
             ELSE
               DVDUMMY(J) = 0.0
             END IF
!                DVDUMMY(J) = DLLJ1(J)
               DVDUMMY(J) = DVDUMMY(J) + 24.0D0 * ( 2.0D0*DLJ1(1,J)*(RHO16*RHO1SQ*RHO1SQ*RHO1-CLJ14(1)/(RHO1SQ*RHO1))&
                   &+14.0D0*DCLJ1(1,J)*(CLJ13(1)/RHO1SQ-CLJ11(1))) !REPULSIVE DERIVATIVE

               DVDUMMY(J) = DVDUMMY(J) + 24.0D0 * (-1.0D0*DLJ1(2,J)*(RHO2SQ*RHO2SQ*RHO2- CLJ8(2)/(RHO2SQ*RHO2))&
                   &- 4.0D0*DCLJ1(2,J)*( CLJ7(2)/RHO2SQ- CLJ5(2))) !ATTRACTIVE DERIVATIVE
!               WRITE(*,*) 'DLJ1(2,J)', DFP(2,J-6),J-6
!               DVDUMMY(J) = DCLJ1(1,J)
            END DO
            DO J=1,3
               TIJ(J) = DVDUMMY(6+J)
               TJI(J) = DVDUMMY(9+J)
            END DO
!        ELSE
                 
!        DO J=1,12
!          DVDUMMY(J)=DRLJ(J)
!        END DO

          ELSE  !NO CUTOFF
             FIJ = 0.0D0
             TIJ = 0.0D0
             TJI = 0.0D0
             DVDUMMY(:)=0.0D0

           IF(LJSITE) THEN
            DO K=1,MAXINTERACTIONS
             DO J=1,3
               DVDUMMY(J) = DVDUMMY(J) + 4.0D0*EPSILON1(K,J1,J2)*(12.0D0*LLJ(11,K)*DLLJ1(K,J))*TERM2(K) + &
                          & 4.0D0*EPSILON1(K,J1,J2)*LLJ(12,K)*TERM3(K)*DRLJ(K,J)! EXTRA LJ SITE DERIVATIVES, CURRENTLY ONLY REPULSIVE
               DVDUMMY(J) = DVDUMMY(J) - ATTR(K)*(4.0D0*EPSILON1(K,J1,J2)*(6.0D0*LLJ(5,K)*DLLJ1(K,J))*TERM2(K) + &
                          & 4.0D0*EPSILON1(K,J1,J2)*LLJ(6,K)*TERM3(K)*DRLJ(K,J)) ! ATTRACTIVE SECONDARY APEX SITE
               FIJ(J) = DVDUMMY(J)
             END DO
             DO J=7,12
               DVDUMMY(J) = DVDUMMY(J) + 4.0D0*EPSILON1(K,J1,J2)*(12.0D0*LLJ(11,K)*DLLJ1(K,J))*TERM2(K) + &
                          & 4.0D0*EPSILON1(K,J1,J2)*LLJ(12,K)*TERM3(K)*DRLJ(K,J)! EXTRA LJ SITE DERIVATIVES, CURRENTLY ONLY REPULSIVE
               DVDUMMY(J) = DVDUMMY(J) - ATTR(K)*(4.0D0*EPSILON1(K,J1,J2)*(6.0D0*LLJ(5,K)*DLLJ1(K,J))*TERM2(K) + &
                          & 4.0D0*EPSILON1(K,J1,J2)*LLJ(6,K)*TERM3(K)*DRLJ(K,J)) ! ATTRACTIVE SECONDARY APEX SITE
             END DO
            END DO
           END IF
!4.0D0*EPSILON0*(12.0D0*LJ11(1)*DLJ1(1,J)-6.0D0*LJ5(2)*DLJ1(2,J))
               FIJ    = FIJ + 24.0D0*(2.D0*RHO112*RHO1*DG1DR - RHO26*RHO2*DG2DR)
               TIJ(1) = DVDF1*DF1PI1 + DVDF2*DF2PI1
               TIJ(2) = DVDF1*DF1PI2 + DVDF2*DF2PI2
               TIJ(3) = DVDF1*DF1PI3 + DVDF2*DF2PI3
               TJI(1) = DVDF1*DF1PJ1 + DVDF2*DF2PJ1
               TJI(2) = DVDF1*DF1PJ2 + DVDF2*DF2PJ2
               TJI(3) = DVDF1*DF1PJ3 + DVDF2*DF2PJ3
               TIJ = DVDUMMY(7:9) + 24.0D0 * TIJ
               TJI = DVDUMMY(10:12) + 24.0D0 * TJI
          END IF

!                G(J3-2)=G(J3-2)+DVDUMMY(1)
!                G(J4-2)=G(J4-2)+DVDUMMY(4)
!                G(J3-1)=G(J3-1)+DVDUMMY(2)
!                G(J4-1)=G(J4-1)+DVDUMMY(5)
!                G(J3  )=G(J3  )+DVDUMMY(3)
!                G(J4  )=G(J4  )+DVDUMMY(6)



!              G(J3-2:J3) = G(J3-2:J3) + DVDUMMY(1:3)
!              G(J4-2:J4) = G(J4-2:J4) + DVDUMMY(4:6)
!              G(J5-2:J5) = G(J5-2:J5) + DVDUMMY(7:9)
!              G(J6-2:J6) = G(J6-2:J6) + DVDUMMY(10:12)

               G(J3-2:J3) = G(J3-2:J3) - FIJ
               G(J4-2:J4) = G(J4-2:J4) + FIJ
               G(J5-2:J5) = G(J5-2:J5) + TIJ
               G(J6-2:J6) = G(J6-2:J6) + TJI

!     END INNER LOOP OVER PARTICLES
124 CONTINUE
            ENDDO

!     END OUTER LOOP OVER PARTICLES

         ENDDO
!        IF(ENERGY<-1.0D5) THEN
!          G=1.0D0
!          ENERGY=0.0D0
!        ELSE
         ENERGY = PYEPSNOT*ENERGY
         G      = PYEPSNOT*G

!        END IF
!998     CONTINUE
        RETURN
      END SUBROUTINE PYGPERIODIC 


SUBROUTINE PARAMONOVNUMFIRSTDER(OLDX,STEST)
USE COMMONS, ONLY : NATOMS
!USE KEY
USE MODHESS
IMPLICIT NONE

DOUBLE PRECISION  :: V(3*NATOMS),ENERGY,X(3*NATOMS),NUMGRAD(3*NATOMS),TRUEGRAD(3*NATOMS),&
        &            OLDX(3*NATOMS),ETEMP(2,3*NATOMS),KSI
INTEGER    :: I,J
LOGICAL    :: GTEST,STEST


KSI=0.00001D0
X(:)=OLDX(:)

!         CALL OLDPARAMONOV(X,TRUEGRAD,ENERGY,GTEST,.FALSE.)
         CALL PYGPERIODIC(X,TRUEGRAD,ENERGY,.TRUE.,.FALSE.)


!WRITE(*,*) 'SF344> IN PARAMONOVSECDER'
DO I=1,3*NATOMS

         X(I)=X(I)-KSI

!         CALL OLDPARAMONOV(X,V,ENERGY,GTEST,.FALSE.)
         CALL PYGPERIODIC(X,V,ENERGY,.TRUE.,.FALSE.)

         ETEMP(1,I)=ENERGY

         X(I)=X(I)+2.0D0*KSI

!         CALL OLDPARAMONOV(X,V,ENERGY,GTEST,.FALSE.)
         CALL PYGPERIODIC(X,V,ENERGY,.TRUE.,.FALSE.)
    
         ETEMP(2,I)=ENERGY
         NUMGRAD(I)=(ETEMP(2,I)-ETEMP(1,I))/(2.0D0*KSI)

         X(I)=X(I)-KSI

WRITE(*,'(3F20.10)') TRUEGRAD(I),NUMGRAD(I), 100.0D0*(TRUEGRAD(I)-NUMGRAD(I))/TRUEGRAD(I)

END DO

STOP
!WRITE(*,*) 'SF344> EXITING AMBERSECDER'
!WRITE(*,*) 'HESSIAN:'
!WRITE(*,'12F10.3') HESS(:,:)
END SUBROUTINE PARAMONOVNUMFIRSTDER

SUBROUTINE PYGSECDER(OLDX,STEST)
USE COMMONS
!USE KEYWORD
USE MODHESS
IMPLICIT NONE

DOUBLE PRECISION  :: V(3*NATOMS),EGB,X(3*NATOMS),OLDX(3*NATOMS),VTEMP(2,3*NATOMS),KSI
INTEGER    :: I,J
LOGICAL    :: GTEST,STEST


KSI=0.00001D0
X(:)=OLDX(:)

!WRITE(*,*) 'SF344> IN PARAMONOVSECDER', NATOMS, X(1)
!KSI=0.0001D0
!X(:)=OLDX(:)

IF (.NOT.ALLOCATED(HESS)) ALLOCATE(HESS(3*NATOMS,3*NATOMS))
HESS(:,:)=0.0D0
!DO I=1,3*NATOMS
!
!         X(I)=X(I)-KSI
!
!         CALL AMBERENERGIES(X,V,ENERGY,GTEST,.FALSE.)
!
!         VTEMP(1,:)=V(:)
!
!         X(I)=X(I)+2.0D0*KSI
!
!         CALL AMBERENERGIES(X,V,ENERGY,GTEST,.FALSE.)
!
!         VTEMP(2,:)=V(:)
!
!                DO J=I,3*NATOM
!                        HESS(I,J)=(VTEMP(2,J)-VTEMP(1,J))/(2.0D0*KSI)
!                        HESS(J,I)=HESS(I,J)
!                END DO
!END DO
V(:)=0.0D0
VTEMP(:,:)=0.0D0

DO I=1,3*NATOMS

!       IF ((I.GT.3*NATOMS/2).AND.(MOD(I,3).EQ.0)) THEN
!        VTEMP(:,:) = 0.0D0
!       ELSE
         X(I)=X(I)-KSI
 
!         CALL OLDPARAMONOV(X,V,EGB,GTEST,.FALSE.)
         CALL PYG(X,V,EGB,.TRUE.,.FALSE.)
!         WRITE(*,*) 'EGB1',EGB
         VTEMP(1,:)=V(:)
 
         X(I)=X(I)+2.0D0*KSI

!         CALL OLDPARAMONOV(X,V,EGB,GTEST,.FALSE.)
         CALL PYG(X,V,EGB,.TRUE.,.FALSE.)
!         WRITE(*,*) 'EGB2',EGB

         VTEMP(2,:)=V(:)

                DO J=I,3*NATOMS
                        HESS(I,J)=(VTEMP(2,J)-VTEMP(1,J))/(2.0D0*KSI)
                        HESS(J,I)=HESS(I,J)
                END DO
!        END IF
END DO
!WRITE(*,*) 'SF344> EXITING PARAMONOVSECDER', NATOMS, X(1)
!WRITE(*,*) 'HESSIAN:'
!WRITE(*,'(12F10.3)') HESS(:,:)

END SUBROUTINE PYGSECDER


SUBROUTINE PYGPERIODICSECDER(OLDX,STEST)
USE COMMONS
!USE KEYWORD
USE MODHESS
IMPLICIT NONE

DOUBLE PRECISION  :: V(3*NATOMS),EGB,X(3*NATOMS),OLDX(3*NATOMS),VTEMP(2,3*NATOMS),KSI
INTEGER    :: I,J
LOGICAL    :: GTEST,STEST


KSI=0.00001D0
X(:)=OLDX(:)

!WRITE(*,*) 'SF344> IN PARAMONOVSECDER', NATOMS, X(1)
!KSI=0.0001D0
!X(:)=OLDX(:)

IF (.NOT.ALLOCATED(HESS)) ALLOCATE(HESS(3*NATOMS,3*NATOMS))
HESS(:,:)=0.0D0
!DO I=1,3*NATOMS
!
!         X(I)=X(I)-KSI
!
!         CALL AMBERENERGIES(X,V,ENERGY,GTEST,.FALSE.)
!
!         VTEMP(1,:)=V(:)
!
!         X(I)=X(I)+2.0D0*KSI
!
!         CALL AMBERENERGIES(X,V,ENERGY,GTEST,.FALSE.)
!
!         VTEMP(2,:)=V(:)
!
!                DO J=I,3*NATOM
!                        HESS(I,J)=(VTEMP(2,J)-VTEMP(1,J))/(2.0D0*KSI)
!                        HESS(J,I)=HESS(I,J)
!                END DO
!END DO
V(:)=0.0D0
VTEMP(:,:)=0.0D0

DO I=1,3*NATOMS

!       IF ((I.GT.3*NATOMS/2).AND.(MOD(I,3).EQ.0)) THEN
!        VTEMP(:,:) = 0.0D0
!       ELSE
         X(I)=X(I)-KSI
 
!         CALL OLDPARAMONOV(X,V,EGB,GTEST,.FALSE.)
         CALL PYGPERIODIC(X,V,EGB,GTEST,.FALSE.)
!         WRITE(*,*) 'EGB1',EGB
         VTEMP(1,:)=V(:)
 
         X(I)=X(I)+2.0D0*KSI

!         CALL OLDPARAMONOV(X,V,EGB,GTEST,.FALSE.)
         CALL PYGPERIODIC(X,V,EGB,GTEST,.FALSE.)
!         WRITE(*,*) 'EGB2',EGB

         VTEMP(2,:)=V(:)

                DO J=I,3*NATOMS
                        HESS(I,J)=(VTEMP(2,J)-VTEMP(1,J))/(2.0D0*KSI)
                        HESS(J,I)=HESS(I,J)
                END DO
!        END IF
END DO
!WRITE(*,*) 'SF344> EXITING PARAMONOVSECDER', NATOMS, X(1)
!WRITE(*,*) 'HESSIAN:'
!WRITE(*,'(12F10.3)') HESS(:,:)

END SUBROUTINE PYGPERIODICSECDER

SUBROUTINE ELLIPSOIDSAATOPOLAR(PX1,PY1,PZ1,ALPHA,BETA,GAMMA,ALPHADEG,BETADEG,GAMMADEG)
! CONVERTS ANGLE-AXIS COORDINATES PX, PY, PZ TO "POLAR-LIKE" ANGLES ALPHA AND BETA
! PX=COS(ALPHA)*COS(ALPHA)
! PY=COS(ALPHA)*SIN(BETA)
! PZ=SIN(ALPHA)
!USE COMMONS
IMPLICIT NONE

DOUBLE PRECISION        :: PX1, PY1, PZ1,PX,PY,PZ

DOUBLE PRECISION,INTENT(OUT)       :: ALPHA, BETA, GAMMA, ALPHADEG, BETADEG, GAMMADEG
DOUBLE PRECISION                   :: PI,TWOPI

! ALPHA: ANGLE OF THE VECTOR WITH THE XY PLANE
! BETA: ANGLE OF THE VECTOR WITH THE X AXIS
! GAMMA: ANGLE OF ROTATION ABOUT THE Z AXIS (ANGLE-AXIS CONVENTION),
!        PROVIDED BY THE MAGNITUDE OF THE VECTOR (SQRT(PX**2+PY**2+PZ**2)-1)
!        SO THAT A VECTOR OF MAGNITUDE 1 MEANS NO ROTATION IS PERFORMED

!IF(PZ*PY*PZ<0) THEN
!       PX = -PX
!       PY = -PY
!       PZ = -PZ
!END IF

PI = ATAN(1.0)*4.0D0
TWOPI = 2.0D0*PI

     GAMMA = SQRT(PX1**2+PY1**2+PZ1**2)

     PX = PX1/GAMMA
     PY = PY1/GAMMA
     PZ = PZ1/GAMMA

     IF(PX.EQ.0.0D0) THEN
       IF(PY>=0.0D0) THEN
        ALPHA = PI/2           ! EULER ANGLE ALPHA 
       ELSE
        ALPHA = -1.0D0*PI/2
       END IF
     ELSE
        IF(PY>=0.0D0) THEN
          IF(PX>0.0D0) THEN
            ALPHA = 1.0D0*ATAN(PY/PX)   ! FIRST QUADRANT
          ELSE    ! PX<0
            ALPHA = 1.0D0*ATAN(PY/PX)+PI       ! SHOULD BE IN THE SECOND QUADRANT
          END IF
        ELSE IF(PY<0.0D0) THEN
          IF(PX>0.0D0) THEN
            ALPHA = 1.0D0*ATAN(PY/PX)    ! FOURTH QUADRANT
          ELSE    ! PX<0
            ALPHA = 1.0D0*ATAN(PY/PX)-PI            ! THIRD QUADRANT
          END IF
        END IF 
     END IF
        BETA = 1.0D0*ACOS(PZ)

ALPHADEG=ALPHA*180/PI
BETADEG=BETA*180/PI
GAMMADEG=GAMMA*180/PI

!WRITE(*,*) 'EXITING ELLIPSOIDSAATOPOLAR'
END SUBROUTINE ELLIPSOIDSAATOPOLAR

      SUBROUTINE AATOEULER(PX,PY,PZ,PHI,THETA,CHI)

      USE COMMONS
      USE KEY

      IMPLICIT NONE
      
      INTEGER          :: I, J, K
      DOUBLE PRECISION :: R(3), P(3), A(3), RM(3,3),PI
      DOUBLE PRECISION, INTENT(OUT) :: PHI, THETA, CHI, PX, PY, PZ
 
      PI = 4.D0 * DATAN(1.D0)

            P(1) = PX
            P(2) = PY
            P(3) = PZ

            CALL ELLIPSOIDROTATION(P, RM)

            PHI   = DATAN2(RM(2,3),RM(1,3)) 
            IF (PHI <= 0.D0) PHI = PHI + 2.D0*PI
            
            THETA = DACOS(RM(3,3))

            CHI   = - DATAN2(RM(3,2),RM(3,1))     
            IF (CHI <= 0.D0) CHI = CHI + 2.D0*PI

            PHI    = PHI*180.D0/PI
            THETA  = THETA*180.D0/PI
            CHI    = CHI*180.D0/PI
           
!            WRITE(3,'(A5,2X,3F20.10,2X,A8,6F20.10)')                 & 
!                 'O', R(1), R(2), R(3),                              &
!                 'ELLIPSE', 2.D0*A(1), 2.D0*A(2), 2.D0*A(3), PHI, THETA, CHI 


      END SUBROUTINE AATOEULER

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE ELLIPSOIDROTATION (P, ROTMAT)
! GIVES BACK THE THREE EULER ANGLES NEEDED FOR DISPLAYING AN ELLIPSOID IN XMAKEMOL
      IMPLICIT NONE

      INTEGER          :: K1, K2
      DOUBLE PRECISION :: E(3,3), I3(3,3), ROTMAT(3,3), P(3), THETA, THETA2

      I3(:,:) = 0.D0
      I3(1,1) = 1.D0
      I3(2,2) = 1.D0
      I3(3,3) = 1.D0

      THETA   = DSQRT(DOT_PRODUCT(P,P))

      IF (THETA == 0.D0) THEN

         ROTMAT = I3

      ELSE

         THETA2  = THETA * THETA
         E(:,:)  = 0.D0
         E(1,2)  = -P(3)
         E(1,3)  =  P(2)
         E(2,3)  = -P(1)
         E(2,1)  = -E(1,2)
         E(3,1)  = -E(1,3)
         E(3,2)  = -E(2,3)
         E       = E/THETA

         ROTMAT = I3 + (1.D0-COS(THETA))*MATMUL(E,E) + E*SIN(THETA)

      ENDIF

      END SUBROUTINE ELLIPSOIDROTATION

SUBROUTINE PYREALIGNXYZ
!USE KEY
USE COMMONS, ONLY : NATOMS
!USE PYMODULE
IMPLICIT NONE

INTEGER NLINES,NFRAMES,J1,J2,J3,NEWNATOMS,REALNATOMS
CHARACTER(LEN=4) ATOMLABELS(NATOMS/2)
CHARACTER(LEN=30) DUMMYLINE
DOUBLE PRECISION COORDSREF(3*NATOMS),COORDSNEW(3*NATOMS),D,DIST2,RMAT(3,3),COORDSDUMMY(3*NATOMS)

 REALNATOMS=NATOMS/2
 OPEN(UNIT=133,FILE='PATH.XYZ',STATUS='OLD') 
 OPEN(UNIT=134,FILE='NEWPATH.XYZ',STATUS='UNKNOWN')

 CALL DETERMINELINES(133,NLINES)
 WRITE(*,*) 'NUMBER OF LINES IN THE PATH.XYZ FILE DETERMINED AS ', NLINES
 NFRAMES=NLINES/(NATOMS/2 + 2)
 WRITE(*,*) 'NUMBER OF FRAMES IN THE PATH.XYZ FILE DETERMINED AS ', NFRAMES

! FIRST FRAME
   READ(133,*) NEWNATOMS
   READ(133,*) DUMMYLINE
   IF(NEWNATOMS/=REALNATOMS) THEN
        WRITE(*,'(A)') 'ERROR - NUMBER OF ATOMS DETERMINED FROM OPTIM DOES NOT MATCH '// &
&                        ' WITH THE NUMBER OF ATOMS DETERMINED FROM THE PATH.XYZ FILE!!!'
        STOP
   END IF
   WRITE(134,*) REALNATOMS
   WRITE(134,*) DUMMYLINE
   DO J2=1,REALNATOMS
        READ(133,'(A1,6X,6G20.10)') ATOMLABELS(J2),COORDSREF(3*J2-2),COORDSREF(3*J2-1),COORDSREF(3*J2), &
        & COORDSREF(3*(J2+REALNATOMS)-2),COORDSREF(3*(J2+REALNATOMS)-1),COORDSREF(3*(J2+REALNATOMS))
        WRITE(134,'(A1,6X,6G20.10)') ATOMLABELS(J2),COORDSREF(3*J2-2),COORDSREF(3*J2-1),COORDSREF(3*J2), &
        & COORDSREF(3*(J2+REALNATOMS)-2),COORDSREF(3*(J2+REALNATOMS)-1),COORDSREF(3*(J2+REALNATOMS))
   END DO


 DO J1=1,NFRAMES-1
   READ(133,*) NEWNATOMS
   READ(133,*) DUMMYLINE
   IF(NEWNATOMS/=REALNATOMS) THEN
        WRITE(*,'(A)') 'ERROR - NUMBER OF ATOMS DETERMINED FROM OPTIM DOES NOT MATCH '// &
&                        ' WITH THE NUMBER OF ATOMS DETERMINED FROM THE PATH.XYZ FILE!!!'
        STOP
   END IF
   DO J2=1,REALNATOMS
        READ(133,'(A1,6X,6G20.10)') ATOMLABELS(J2),COORDSNEW(3*J2-2),COORDSNEW(3*J2-1),COORDSNEW(3*J2), &
        & COORDSNEW(3*(J2+REALNATOMS)-2),COORDSNEW(3*(J2+REALNATOMS)-1),COORDSNEW(3*(J2+REALNATOMS))
   END DO
   COORDSDUMMY(:)=COORDSNEW(:)
!   COORDSDUMMY(50)=COORDSNEW(50)*10.0D0
   CALL MINPERMDIST(COORDSREF,COORDSDUMMY,NATOMS,.FALSE.,1.0,1.0,1.0,.FALSE.,.FALSE.,D,DIST2,.FALSE.,RMAT)
   D=D**2 ! SINCE MINPERMDIST NOW RETURNS THE DISTANCE
   COORDSNEW(:)=COORDSDUMMY(:)
   WRITE(*,*) 'WRITING FRAME ', J1+1
   WRITE(134,*) REALNATOMS
   WRITE(134,*) DUMMYLINE
   DO J2=1,REALNATOMS
     WRITE(134,'(A1,6X,6G20.10)') ATOMLABELS(J2),COORDSNEW(3*J2-2),COORDSNEW(3*J2-1),COORDSNEW(3*J2), &
        & COORDSNEW(3*(J2+REALNATOMS)-2),COORDSNEW(3*(J2+REALNATOMS)-1),COORDSNEW(3*(J2+REALNATOMS))
   END DO
   COORDSREF(:)=COORDSNEW(:)
 END DO


END SUBROUTINE PYREALIGNXYZ

SUBROUTINE PYRANDOMSWAP(COORDSREF,COORDSDUMMY,D)
USE COMMONS, ONLY : NATOMS
USE KEY, ONLY : PYBINARYTYPE1
IMPLICIT NONE

INTEGER :: NP, RANDOM1, RANDOM2, J1, J2, J3
DOUBLE PRECISION :: COORDSSTORE(3,4),DPRAND,COORDSREF(3*NATOMS),COORDSDUMMY(3*NATOMS),D
DOUBLE PRECISION :: RMAT(3,3),BESTCOORDS(3*NATOMS),DBEST

DBEST=1.0D10
BESTCOORDS(:)=COORDSDUMMY(:)
DO J3=1,50

DO J2=1,1000000
COORDSDUMMY(:)=BESTCOORDS(:)
! FIRST SWAP 
        DO
          RANDOM1=INT(DPRAND()*100)
          IF(RANDOM1<=PYBINARYTYPE1.AND.RANDOM1>=1) EXIT
        END DO
        DO
          RANDOM2=INT(DPRAND()*100)
          IF(RANDOM2<=PYBINARYTYPE1.AND.RANDOM2>=1) EXIT
        END DO

        !RANDOM1=PYSWAP(1)
        !RANDOM2=PYSWAP(2)
        COORDSSTORE(:,:)=0.0D0
        ! SELECT A BODY FROM THE TWO TYPES
        DO J1=1,3
         COORDSSTORE(J1,1)=COORDSDUMMY(3*(RANDOM1-1)+J1)
         COORDSSTORE(J1,2)=COORDSDUMMY(3*(RANDOM2-1)+J1)
         COORDSSTORE(J1,3)=COORDSDUMMY(3*NATOMS/2+3*(RANDOM1-1)+J1)
         COORDSSTORE(J1,4)=COORDSDUMMY(3*NATOMS/2+3*(RANDOM2-1)+J1)
        END DO

        ! SWAP COORDINATES
        !WRITE(MYUNIT,*) COORDSSTORE(:,1)
        !WRITE(MYUNIT,*) COORDSSTORE(:,3)
        !WRITE(MYUNIT,*) COORDSSTORE(:,2)
        !WRITE(MYUNIT,*) COORDSSTORE(:,4)
        DO J1=1,3
         COORDSDUMMY(3*(RANDOM1-1)+J1)=COORDSSTORE(J1,2)
         COORDSDUMMY(3*(RANDOM2-1)+J1)=COORDSSTORE(J1,1)
         COORDSDUMMY(3*NATOMS/2+3*(RANDOM1-1)+J1)=COORDSSTORE(J1,4)
         COORDSDUMMY(3*NATOMS/2+3*(RANDOM2-1)+J1)=COORDSSTORE(J1,3)
        END DO

        DO
          RANDOM1=INT(DPRAND()*100)
          IF(RANDOM1>PYBINARYTYPE1.AND.RANDOM1<=NATOMS/2) EXIT
        END DO
        DO
          RANDOM2=INT(DPRAND()*100)
          IF(RANDOM2>PYBINARYTYPE1.AND.RANDOM2<=NATOMS/2) EXIT
        END DO

        !RANDOM1=PYSWAP(1)
        !RANDOM2=PYSWAP(2)
        COORDSSTORE(:,:)=0.0D0
        ! SELECT A BODY FROM THE TWO TYPES
        DO J1=1,3
         COORDSSTORE(J1,1)=COORDSDUMMY(3*(RANDOM1-1)+J1)
         COORDSSTORE(J1,2)=COORDSDUMMY(3*(RANDOM2-1)+J1)
         COORDSSTORE(J1,3)=COORDSDUMMY(3*NATOMS/2+3*(RANDOM1-1)+J1)
         COORDSSTORE(J1,4)=COORDSDUMMY(3*NATOMS/2+3*(RANDOM2-1)+J1)
        END DO

        ! SWAP COORDINATES
        !WRITE(MYUNIT,*) COORDSSTORE(:,1)
        !WRITE(MYUNIT,*) COORDSSTORE(:,3)
        !WRITE(MYUNIT,*) COORDSSTORE(:,2)
        !WRITE(MYUNIT,*) COORDSSTORE(:,4)
        DO J1=1,3
         COORDSDUMMY(3*(RANDOM1-1)+J1)=COORDSSTORE(J1,2)
         COORDSDUMMY(3*(RANDOM2-1)+J1)=COORDSSTORE(J1,1)
         COORDSDUMMY(3*NATOMS/2+3*(RANDOM1-1)+J1)=COORDSSTORE(J1,4)
         COORDSDUMMY(3*NATOMS/2+3*(RANDOM2-1)+J1)=COORDSSTORE(J1,3)
        END DO

   IF(MOD(J2,1000)==0)   THEN 
   CALL MINPERMDIST(COORDSREF,COORDSDUMMY,NATOMS,.FALSE.,1.0,1.0,1.0,.FALSE.,.FALSE.,D,D,.FALSE.,RMAT)
   D=D**2
        IF(D<DBEST) THEN
          DBEST=D
          BESTCOORDS(:)=COORDSDUMMY(:)
!        ELSE IF (D-DBEST<0.1D0.AND.DPRAND()<0.5D0) THEN
!          DBEST=D
!          BESTCOORDS(:)=COORDSDUMMY(:)
        END IF
   END IF
END DO
WRITE(*,*) 'AFTER 10000 SWAPS, DBEST= ', DBEST

END DO

END SUBROUTINE PYRANDOMSWAP

SUBROUTINE DETERMINELINES(NUNIT,NLINES)
IMPLICIT NONE
INTEGER NUNIT, NLINES, IOSTATUS
CHARACTER(LEN=10) CHECK

REWIND(NUNIT)

NLINES=0
DO
  IF(IOSTATUS<0) EXIT
  NLINES = NLINES + 1
  READ(NUNIT,*,IOSTAT=IOSTATUS) CHECK
!  WRITE(*,*) CHECK,NUNIT
END DO
  NLINES = NLINES - 1
  REWIND(NUNIT)
RETURN


END SUBROUTINE DETERMINELINES

