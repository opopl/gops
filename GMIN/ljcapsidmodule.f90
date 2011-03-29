
MODULE LJCAPSIDMODULE

 INTEGER          :: I, J, J1, J2, J3, J4, J5, J6, K1, K2, OFFSET, REALNATOMS
 DOUBLE PRECISION, ALLOCATABLE :: RMIVEC(:,:,:), DPI1RMVEC(:,:,:), DPI2RMVEC(:,:,:), DPI3RMVEC(:,:,:)
 DOUBLE PRECISION, ALLOCATABLE :: PSCALEFAC1VEC(:),PSCALEFAC2VEC(:),EPSILON1(:,:,:),AEZR1(:,:,:), AEZR2(:,:,:)

 DOUBLE PRECISION :: ANGLE,ANGLE2,PI,TWOPI,SIGMA1,CUT,RON,RON2,RANGE2INV3

 DOUBLE PRECISION ::I3(3,3) 

END MODULE LJCAPSIDMODULE

SUBROUTINE INITIALISEPYGPERIODIC

USE COMMONS, ONLY: BOXLX,BOXLY,BOXLZ,PARAMONOVPBCX,PARAMONOVPBCY,PARAMONOVPBCZ,PCUTOFF,PARAMONOVCUTOFF,&
                &       NATOMS,PYA1BIN,PYA2BIN,PYSIGNOT,PYEPSNOT,RADIFT,LJSITE,BLJSITE,PEPSILON1,&
                &       PSCALEFAC1,PSCALEFAC2,MAXINTERACTIONS,PYBINARYT,PYBINARYTYPE1,MYUNIT,VT, &
                &       PEPSILONATTR, PSIGMAATTR, LJSITEATTR, LJSITECOORDST, LJSITECOORDS

USE PYMODULE

IMPLICIT NONE
   WRITE(MYUNIT,*) 'INITIALISING VARIABLES FOR PY',NATOMS 
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
        WRITE(MYUNIT,'(A,3F8.3)') 'REPULSIVE LJ SITE COORDINATES WILL BE ', LJSITECOORDS(:)
    END IF
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

END SUBROUTINE INITIALISEPYGPERIODIC

! PY POTENTIAL, DC430'S IMPLEMENTATION
! WITH PBC AND CONTINUOUS CUTOFF ADDED
! PLUS EXTRA LJ SITE
      SUBROUTINE PYGPERIODIC (X, G, ENERGY, GTEST)

       USE COMMONS, ONLY: BOXLX,BOXLY,BOXLZ,PARAMONOVPBCX,PARAMONOVPBCY,PARAMONOVPBCZ,PCUTOFF,PARAMONOVCUTOFF,&
                &       NATOMS,PYA1BIN,PYA2BIN,PYSIGNOT,PYEPSNOT,RADIFT,LJSITE,BLJSITE,PEPSILON1,& 
                &       PSCALEFAC1,PSCALEFAC2,MAXINTERACTIONS,PYBINARYT,PYBINARYTYPE1,MYUNIT,VT,PEPSILONATTR,PSIGMAATTR,&
                &       FROZEN

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

      DO J1=1,REALNATOMS
        J2=3*J1

        IF (PARAMONOVPBCX) THEN
!ENSURE X COMPONENT OF PARTICLE 1 VECTOR IS WITHIN BOXLX/2 OF ZERO. IF IT ISN'T THEN SUBTRACT INTEGER NUMBER OF BOXLX'S SUCH THAT IT IS.
                X(J2-2)=X(J2-2)-BOXLX*NINT(X(J2-2)/BOXLX)
        ENDIF

        IF (PARAMONOVPBCY) THEN
!ENSURE Y COMPONENT OF PARTICLE 1 VECTOR IS WITHIN BOXLY/2 OF ZERO. IF IT ISN'T THEN SUBTRACT INTEGER NUMBER OF BOXLY'S SUCH THAT IT IS.
                X(J2-1)=X(J2-1)-BOXLY*NINT(X(J2-1)/BOXLY)
        END IF

        IF (PARAMONOVPBCZ) THEN
!ENSURE Z COMPONENT OF PARTICLE 1 VECTOR IS WITHIN BOXLZ/2 OF ZERO. IF IT ISN'T THEN SUBTRACT INTEGER NUMBER OF BOXLZ'S SUCH THAT IT IS.
                X(J2)=X(J2)-BOXLZ*NINT(X(J2)/BOXLZ)
        ENDIF
         
      END DO

         ENERGY = 0.D0
         G(:)   = 0.D0
        
        DO J1=1, REALNATOMS
            J3      = 3*J1
            J5      = OFFSET + J3
            RI      = X(J3-2:J3)
            P       = X(J5-2:J5)
            ANGLE=SQRT(DOT_PRODUCT(P,P))
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

            CALL RMDRVT(P, RMIVEC(J1,:,:), DPI1RMVEC(J1,:,:), DPI2RMVEC(J1,:,:), DPI3RMVEC(J1,:,:), GTEST)

        END DO

         DO J1 = 1, REALNATOMS
          
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
                 
!        IF(FROZEN(J1)) THEN 
                VDUMMY=0.0D0
!        ELSE
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
                 VDUMMY=VDUMMY+4.0D0*EPSILON1(K,J1,J2)*TERM2(K)*(LLJ(12,K) - ATTR(K)*LLJ(6,K)) ! EXTRA LJ SITES (12-6)
!                 VDUMMY=VDUMMY+4.0D0*EPSILON1(K,J1,J2)*TERM2(K)*LLJ(6,K) ! EXTRA LJ SITES (REPULSIVE LJ6, TESTING)
              END DO ! K=1,4
             END IF ! IF(LJSITE)

        IF(PARAMONOVCUTOFF) THEN
                !REPULSIVE POTENTIAL AND PERIODIC CUTOFF CORRECTIONS
                VDUMMY = VDUMMY + 4.0D0 * PYEPSNOT * ( RHO112 + 6.0D0*CLJ12(1)*CLJ2(1)/RHO1SQ-7.0D0*CLJ12(1))
                    !                      LJ12(1)+(6.0D0*CLJ12(1))*CLJ2(1)/LJ2(1)-7.0D0*CLJ12(1)
                !ATTRACTIVE POTENTIAL AND PERIODIC CUTOFF CORRECTIONS
                VDUMMY = VDUMMY + 4.0D0 * PYEPSNOT * (- RHO26 - 3.0D0* CLJ6(2)*CLJ2(2)/RHO2SQ+4.0D0* CLJ6(2))
                    !                      (-LJ6(2)+(-3.0D0*CLJ6(2))*CLJ2(2)/LJ2(2) +4.0D0*CLJ6(2))
        ELSE
                VDUMMY = VDUMMY + 4.0D0 * PYEPSNOT * (RHO112 - 1.0D0 * RHO26)
        END IF

!       END IF ! IF(FROZEN)

               ENERGY = ENERGY + VDUMMY
               VT(J1) = VT(J1) + VDUMMY
               VT(J2) = VT(J2) + VDUMMY        ! PAIR POTENTIALS


            DVDUMMY(:) = 0.0D0

!     CALCULATE GRADIENT
             FIJ = 0.0D0
             TIJ = 0.0D0
             TJI = 0.0D0

 !       IF(.NOT.FROZEN(J1)) THEN
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
               DVDUMMY(J) = DVDUMMY(J) + 24.0D0 * PYEPSNOT * ( &
               & 2.0D0*RHO1SQ*DG1DR(J)*(RHO16*RHO1SQ*RHO1SQ*RHO1-CLJ14(1)/(RHO1SQ*RHO1))&
               &+14.0D0*DCLJ1(1,J)*(CLJ13(1)/RHO1SQ-CLJ11(1))) 
               !REPULSIVE DERIVATIVE

               DVDUMMY(J) = DVDUMMY(J) + 24.0D0 * PYEPSNOT * (-1.0D0*RHO2SQ*DG2DR(J)*(RHO2SQ*RHO2SQ*RHO2-CLJ8(2)/(RHO2SQ*RHO2))&
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
               DVDUMMY(J) = DVDUMMY(J) + 24.0D0 * PYEPSNOT * ( 2.0D0*DLJ1(1,J)*(RHO16*RHO1SQ*RHO1SQ*RHO1-CLJ14(1)/(RHO1SQ*RHO1))&
                   &+14.0D0*DCLJ1(1,J)*(CLJ13(1)/RHO1SQ-CLJ11(1))) !REPULSIVE DERIVATIVE

               DVDUMMY(J) = DVDUMMY(J) + 24.0D0 * PYEPSNOT * (-1.0D0*DLJ1(2,J)*(RHO2SQ*RHO2SQ*RHO2- CLJ8(2)/(RHO2SQ*RHO2))&
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
                          & 4.0D0*EPSILON1(K,J1,J2)*LLJ(12,K)*TERM3(K)*DRLJ(K,J)! EXTRA LJ SITE DERIVATIVES, CURRENTLY ONLY REPULSIVE (LJ12)
!               DVDUMMY(J) = DVDUMMY(J) + 4.0D0*EPSILON1(K,J1,J2)*(6.0D0*LLJ(5,K)*DLLJ1(K,J))*TERM2(K)! + &
!!                          & 4.0D0*EPSILON1(K,J1,J2)*LLJ(12,K)*TERM3(K)*DRLJ(K,J)! EXTRA LJ SITE DERIVATIVES, CURRENTLY ONLY REPULSIVE (LJ6)
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
           END IF !IF LJSITE 
!        END IF !IF(.NOT.FROZEN)
!4.0D0*EPSILON0*(12.0D0*LJ11(1)*DLJ1(1,J)-6.0D0*LJ5(2)*DLJ1(2,J))
               FIJ    = FIJ + 24.0D0*PYEPSNOT*(2.D0*RHO112*RHO1*DG1DR - 1.0D0 * RHO26*RHO2*DG2DR)
               TIJ(1) = DVDF1*DF1PI1 + 1.0D0* DVDF2*DF2PI1
               TIJ(2) = DVDF1*DF1PI2 + 1.0D0*  DVDF2*DF2PI2
               TIJ(3) = DVDF1*DF1PI3 + 1.0D0*  DVDF2*DF2PI3
               TJI(1) = DVDF1*DF1PJ1 + 1.0D0*  DVDF2*DF2PJ1
               TJI(2) = DVDF1*DF1PJ2 + 1.0D0*  DVDF2*DF2PJ2
               TJI(3) = DVDF1*DF1PJ3 + 1.0D0*  DVDF2*DF2PJ3
               TIJ = DVDUMMY(7:9) + 24.0D0 * PYEPSNOT * TIJ
               TJI = DVDUMMY(10:12) + 24.0D0 * PYEPSNOT * TJI
          END IF ! CUTOFF/NO CUTOFF BLOCK

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
!         ENERGY = PYEPSNOT*ENERGY
!         G      = PYEPSNOT*G
!        END IF
!998     CONTINUE
        DO J1=1,NATOMS
                J2=3*J1
                IF(FROZEN(J1)) THEN
                        G(J2-2)=0.0D0
                        G(J2-1)=0.0D0
                        G(J2)=0.0D0
                END IF
        END DO
        RETURN
      END SUBROUTINE PYGPERIODIC 



! PY POTENTIAL, SF344'S IMPLEMENTATION
       SUBROUTINE OLDPARAMONOV(X,V,EGB,GTEST,STEST)
       USE COMMONS, ONLY: BOXLX,BOXLY,BOXLZ,PARAMONOVPBCX,PARAMONOVPBCY,PARAMONOVPBCZ,PCUTOFF,PARAMONOVCUTOFF,&
                &       NATOMS,PSIGMA0,PEPSILON0,PARAMA1,PARAMB1,PARAMC1,PARAMA2,PARAMB2,PARAMC2,VT,LJSITE,PEPSILON1
       !USE KEYWORD
       IMPLICIT NONE
!      ANISOTROPIC POTENTIAL BY PARAMONOV ET. AL., JCP 123, 194111 (2005)
!      FINALLY WORKING FOR TRIAXIAL SYSTEMS AS WELL

       DOUBLE PRECISION :: X1,Y1,Z1,X2,Y2,Z2,A1(2),B1(2),C1(2),F1(2),DF1(2,12),R,RLJ,RLJ2,XLJ(2,3),DRLJ(12), &
                        &  TERM2, TERM3,RON2,RON,CUT,RANGE2INV3
       DOUBLE PRECISION :: PX1(2),PY1(2),PZ1(2),ALPHA(2),BETA(2),GAMMA(2), &
                        &                   SA(2),SB(2),SG(2),CA(2),CB(2),CG(2),&
                        &                   PX(2),PY(2),PZ(2),QX(2),QY(2),QZ(2),RX(2),RY(2),RZ(2)
       DOUBLE PRECISION :: DA(2,3),DB(2,3),DG(2,3),DPX(2,3),DPY(2,3),DPZ(2,3),DQX(2,3),DQY(2,3),&
                        &                   DQZ(2,3),DRX(2,3),DRY(2,3),DRZ(2,3)
       DOUBLE PRECISION :: LJ1(2),LJ2(2),LJ3(2),LJ5(2),LJ6(2), &
                        &                  LJ12(2), LJ4(2), LJ7(2), LJ11(2),DLJ1(2,12),DR(6)
       DOUBLE PRECISION :: LAMBDA(2), A11(3,3),B11(3,3),C_1(6),C1INV(3,3),W(3),DUMMY,STEPSIZE,DF_LAMBDAOLD,LAMBDAOLD
       DOUBLE PRECISION :: DUMMY1, DUMMY2, S_DUMMY1, C_DUMMY1, S_DUMMY2, C_DUMMY2, VDUMMY, DVDUMMY(12),EPSILON0
       DOUBLE PRECISION :: DUMMYX, DUMMYY,DUMMYZ, SIGMA0(2), DF_LAMBDA(1), EGB,X(3*NATOMS), &
                        &  V(3*NATOMS),DPRAND,CONVG,U1(3),U2(3),R1(3),R1T(1,3)
       INTEGER   :: I,J,K,L,ITER=0,REALNATOMS,J1,J2,J3,K2,J4
       LOGICAL   :: GTEST,STEST
       DOUBLE PRECISION :: G_INVR(3),G_INVRT(1,3),S_LAMBDA(1), RANDOM,PI,FCHECK(3),DCHECK
       DOUBLE PRECISION :: DGU1(3,3,3),DGU2(3,3,3),DXU1(3,3),DXU2(3,3),DFU1(3,1),DFU2(3,1)
       DOUBLE PRECISION :: DAA(3,2,3,3),DBB(3,2,3,3)

       DOUBLE PRECISION :: LLJ1,LLJ2,LLJ3,LLJ4,LLJ5,LLJ6,LLJ7,LLJ11,LLJ12,DLLJ1(12),SIGMA1,EPSILON1
       DOUBLE PRECISION :: CLJ1(2),CLJ2(2),CLJ3(2),CLJ4(2),CLJ5(2),CLJ6(2)
       DOUBLE PRECISION :: CLJ7(2),CLJ8(2),CLJ11(2),CLJ12(2),CLJ13(2),CLJ14(2),DCLJ1(2,12)

       EPSILON1=PEPSILON1(1)   ! WE ARE GOING TO USE SIGMA1 AND EPSILON1 FOR THE EXTRA LJ SITES
       SIGMA1=1.0D0!PSIGMA1 IS NONEXISTENT FROM NOW ON
       REALNATOMS=NATOMS/2
       PI=4.0D0*ATAN(1.0)
! CUTOFF STUFF FOR THE EXTRA LJ SITES
       CUT = PCUTOFF*PCUTOFF ! CUTOFF SQUARED
       RON = PCUTOFF*0.9D0
       RON2 = RON*RON
       RANGE2INV3=1.0D0/(CUT-RON2)**3
       TERM2=1.0D0
       TERM3=0.0D0
!      VT(1:REALNATOMS)=0.0D0
       V(1:6*REALNATOMS)=0.0D0
!      WRITE(*,*) "IN ENERGY CYCLE"
!      WRITE(*,*) X(:)
!      WRITE(*,*) BOXLX, BOXLY, BOXLZ, PCUTOFF
!      WRITE (*,*) PARAMONOVPBCX,PARAMONOVPBCY,PARAMONOVPBCZ, PARAMONOVCUTOFF

       EGB=0.0D0

       SIGMA0(1)=PSIGMA0(1)
       SIGMA0(2)=PSIGMA0(2)
       EPSILON0=PEPSILON0

       A1(1)=PARAMA1**2
       B1(1)=PARAMB1**2
       C1(1)=PARAMC1**2

       A1(2)=PARAMA2**2
       B1(2)=PARAMB2**2
       C1(2)=PARAMC2**2

!      WRITE(*,*) A1(:),B1(:),C1(:)

       DO J1=1,REALNATOMS
        J2=3*J1

        IF (PARAMONOVPBCX) THEN
                !ENSURE X COMPONENT OF PARTICLE 1 VECTOR IS WITHIN BOXLX/2 OF ZERO. IF IT ISN'T THEN SUBTRACT INTEGER NUMBER OF BOXLX'S SUCH THAT IT IS.
                X(J2-2)=X(J2-2)-BOXLX*NINT(X(J2-2)/BOXLX)
        ENDIF

        IF (PARAMONOVPBCY) THEN
                !ENSURE Y COMPONENT OF PARTICLE 1 VECTOR IS WITHIN BOXLY/2 OF ZERO. IF IT ISN'T THEN SUBTRACT INTEGER NUMBER OF BOXLY'S SUCH THAT IT IS.
                X(J2-1)=X(J2-1)-BOXLY*NINT(X(J2-1)/BOXLY)
        END IF

        IF (PARAMONOVPBCZ) THEN
                !ENSURE Z COMPONENT OF PARTICLE 1 VECTOR IS WITHIN BOXLZ/2 OF ZERO. IF IT ISN'T THEN SUBTRACT INTEGER NUMBER OF BOXLZ'S SUCH THAT IT IS.
                X(J2)=X(J2)-BOXLZ*NINT(X(J2)/BOXLZ)
        ENDIF


         
       END DO

       DO J1=1,REALNATOMS
        J2=3*J1
        
        X1 = X(J2-2)
        Y1 = X(J2-1)
        Z1 = X(J2  )

!        PX, PY, PZ: ANGLE-AXIS COORDINATE SYSTEM (BUT: RIGHT HAND ROTATION!)
        PX1(1) = X(3*REALNATOMS+J2-2)
        PY1(1) = X(3*REALNATOMS+J2-1)
        PZ1(1) = X(3*REALNATOMS+J2  )
!        IF(PX1(1)==0.0D0) PX1(1)=PX1(1)+4.0D0*PI**2
!        IF(PY1(1)==0.0D0) PY1(1)=PY1(1)+4.0D0*PI**2
!        IF(PZ1(1)==0.0D0) PZ1(1)=PZ1(1)+4.0D0*PI**2

        DO J3=J1+1,REALNATOMS
                VDUMMY=0.0D0
                K2=3*J3
!                IF (PARAMONOVPBCX) THEN
!                   !ENSURE X COMPONENT OF PARTICLE 2 VECTOR IS WITHIN BOXLX/2 OF ZERO. IF IT ISN'T THEN SUBTRACT INTEGER NUMBER OF BOXLX'S SUCH THAT IT IS.
!                        X(K2-2)=X(K2-2)-BOXLX*NINT(X(K2-2)/BOXLX)
!                ENDIF
!
!                IF (PARAMONOVPBCY) THEN
!                   !ENSURE Y COMPONENT OF PARTICLE 2 VECTOR IS WITHIN BOXLY/2 OF ZERO. IF IT ISN'T THEN SUBTRACT INTEGER NUMBER OF BOXLY'S SUCH THAT IT IS.
!                        X(K2-1)=X(K2-1)-BOXLY*NINT(X(K2-1)/BOXLY)
!                ENDIF
!
!                IF (PARAMONOVPBCZ) THEN
!                   !ENSURE Z COMPONENT OF PARTICLE 2 VECTOR IS WITHIN BOXLZ/2 OF ZERO. IF IT ISN'T THEN SUBTRACT INTEGER NUMBER OF BOXLZ'S SUCH THAT IT IS.
!                        X(K2)=X(K2)-BOXLZ*NINT(X(K2)/BOXLZ)
!                ENDIF
                X2 = X(K2-2)
                Y2 = X(K2-1)
                Z2 = X(K2  )
                PX1(2) = X(3*REALNATOMS+K2-2)
                PY1(2) = X(3*REALNATOMS+K2-1)
                PZ1(2) = X(3*REALNATOMS+K2  )
!                IF(PX1(2)==0.0D0) PX1(2)=PX1(2)+4.0D0*PI**2
!                IF(PY1(2)==0.0D0) PY1(2)=PY1(2)+4.0D0*PI**2
!                IF(PZ1(2)==0.0D0) PZ1(2)=PZ1(2)+4.0D0*PI**2

!            DERIVATIVE CHECK
!             DO L=1,3
 
                R1(1)=X2-X1
                R1(2)=Y2-Y1
                R1(3)=Z2-Z1

                !COMPUTE VECTOR TO IMAGE OF SECOND PARTICLE THAT IS WITHIN BOXLX/2 OF FIRST PARTICLE
                IF (PARAMONOVPBCX) THEN
                        !WRITE(*,*) "X1,X2,R1(1)"
                        !WRITE(*,*) X1,X2,R1(1)
                        !FIND THE X-COMPONENT OF THE VECTOR TO THE NEAREST IMAGE OF THE SECOND PARTICLE
                        R1(1)=R1(1)-BOXLX*NINT(R1(1)/BOXLX)
                        
                        !RECOMPUTE X2 WITH THE NEW COMPONENT. DON'T WORRY IF IT'S OUTSIDE THE BOX.
                        X2=R1(1)+X1
                        !WRITE(*,*) "X1,X2,R1(1)"
                        !WRITE(*,*) X1,X2,R1(1)
                ENDIF

                !SET PERIODIC BOUNDARY CONDITIONS FOR Y DIMENSION.
                IF (PARAMONOVPBCY) THEN
                        !WRITE(*,*) "Y1,Y2,R1(2)"
                        !WRITE(*,*) Y1,Y2,R1(2)

                        !FIND THE Y-COMPONENT OF THE VECTOR TO THE NEAREST IMAGE OF THE SECOND PARTICLE
                        R1(2)=R1(2)-BOXLY*NINT(R1(2)/BOXLY)

                        !RECOMPUTE Y2 WITH THE NEW COMPONENT. DON'T WORRY IF IT'S OUTSIDE THE BOX.
                        Y2=R1(2)+Y1
                        !WRITE(*,*) "Y1,Y2,R1(2)"
                        !WRITE(*,*) Y1,Y2,R1(2)
                ENDIF

                !SET PERIODIC BOUNDARY CONDITIONS FOR THE Z DIMENSION
                IF (PARAMONOVPBCZ) THEN
                        !WRITE(*,*) "Z1,Z2,R1(3)"
                        !WRITE(*,*) Z1,Z2,R1(3)
                        !FIND THE Z-COMPONENT OF THE VECTOR TO THE NEAREST IMAGE OF THE SECOND PARTICLE
                        R1(3)=R1(3)-BOXLZ*NINT(R1(3)/BOXLZ)
                        !RECOMPUTE Z2 WITH THE NEW COMPONENT.
                        Z2=R1(3)+Z1
                        !WRITE(*,*) "Z1,Z2,R1(3)"
                        !WRITE(*,*) Z1,Z2,R1(3)
                END IF

                !CHECK FOR COLD FUSION
                IF(SQRT((X2-X1)**2+(Y2-Y1)**2+(Z2-Z1)**2).LE.0.01D0) THEN
                        WRITE(*,*) 'SF344> COLD FUSION DETECTED'
                        VDUMMY=-1.0D20
                        EGB=-1.0D20
                        V(:)=0.0D0
                        GOTO 999
                END IF
                                
                R=SQRT((X2-X1)**2+(Y2-Y1)**2+(Z2-Z1)**2)

!IF WE ARE CUTTING OFF THE POTENTIAL AND THE INTER PARTICLE DISTANCE IS >=PCUTOFF DON'T BOTHER DOING ANY MATHS. THE ANSWER IS ZERO.
                IF ((.NOT.PARAMONOVCUTOFF).OR.(PARAMONOVCUTOFF.AND.(R<PCUTOFF))) THEN
                        R1T(1,1)=R1(1)
                        R1T(1,2)=R1(2)
                        R1T(1,3)=R1(3)

                        DR(1)=-1.0D0/R*(X2-X1)
                        DR(2)=-1.0D0/R*(Y2-Y1)
                        DR(3)=-1.0D0/R*(Z2-Z1)
                        DR(4)=-DR(1)
                        DR(5)=-DR(2)
                        DR(6)=-DR(3)
                          
                        DO K=1,2

    
                                STEPSIZE=1.0D0
                                LAMBDA(K)=0.5D0
                                IF (K==2) THEN
                                        LAMBDA(2)=LAMBDA(1)
                                END IF
                                ITER=0
                                CONVG=1.0D-8
                                DF_LAMBDAOLD=1.0D6
 
                                DO I=1,2
                                        GAMMA(I) = SQRT(PX1(I)**2+PY1(I)**2+PZ1(I)**2)
!                                        WRITE(*,*) 'GAMMA(I) BEFORE',I,GAMMA(I)
!                                        GAMMA(I) = GAMMA(I)-2.0D0*PI*NINT(GAMMA(I)/(2.0D0*PI))+2.0D0*PI
!                                        WRITE(*,*) 'GAMMA(I) AFTER',I,GAMMA(I)

                                        DG(I,1) = PX1(I)/GAMMA(I)
                                        DG(I,2) = PY1(I)/GAMMA(I)
                                        DG(I,3) = PZ1(I)/GAMMA(I)
                                        PX(I) = PX1(I)/GAMMA(I)
                                        PY(I) = PY1(I)/GAMMA(I)
                                        PZ(I) = PZ1(I)/GAMMA(I)
ALPHA(I)=ATAN2(PY(I),PX(I))
!                                                                DA(I,1) = -1.0D0*PY1(I)/(PX1(I)**2+PY1(I)**2)
!                                                                DA(I,2) = PX1(I)/(PX1(I)**2+PY1(I)**2)

                                                                DA(I,1) = -1.0D0*PY1(I)/(PX1(I)**2+PY1(I)**2)
                                                                DA(I,2) = PX1(I)/(GAMMA(I)**2-PZ1(I)**2)
                                                                DA(I,3) = 0.0D0
 
!IF(0) THEN
                                        IF(PX(I).EQ.0.0D0) THEN
                                                IF(PY(I)>=0.0D0) THEN
                                                        ALPHA(I) = PI/2           ! EULER ANGLE ALPHA 
                                                ELSE
                                                        ALPHA(I) = -1.0D0*PI/2
                                                END IF
                                        ELSE
                                                IF(PY(I)>=0.0D0) THEN
                                                        IF(PX(I)>0.0D0) THEN
                                                                ALPHA(I) = 1.0D0*ATAN(PY(I)/PX(I))   ! FIRST QUADRANT
                                                                DA(I,1) = -1.0D0*PY1(I)/(PX1(I)**2+PY1(I)**2)
                                                                DA(I,2) = PX1(I)/(PX1(I)**2+PY1(I)**2)
                                                                DA(I,3) = 0.0D0
                                                        ELSE    ! PX<0
                                                                ALPHA(I) = 1.0D0*ATAN(PY(I)/PX(I))+PI ! SHOULD BE IN THE SECOND QUADRANT
                                                                DA(I,1) = -1.0D0*PY1(I)/(PX1(I)**2+PY1(I)**2)
                                                                DA(I,2) = PX1(I)/(PX1(I)**2+PY1(I)**2)
                                                                DA(I,3) = 0.0D0
                                                        END IF
                                                ELSE IF(PY(I)<0.0D0) THEN
                                                        IF(PX(I)>0.0D0) THEN
                                                                ALPHA(I) = 1.0D0*ATAN(PY(I)/PX(I))    ! FOURTH QUADRANT
                                                                DA(I,1) = -1.0D0*PY1(I)/(PX1(I)**2+PY1(I)**2)
                                                                DA(I,2) = PX1(I)/(PX1(I)**2+PY1(I)**2)
                                                                DA(I,3) = 0.0D0
                                                        ELSE    ! PX<0
                                                                ALPHA(I) = 1.0D0*ATAN(PY(I)/PX(I))-PI            ! THIRD QUADRANT
                                                                DA(I,1) = -1.0D0*PY1(I)/(PX1(I)**2+PY1(I)**2)
                                                                DA(I,2) = PX1(I)/(PX1(I)**2+PY1(I)**2)
                                                                DA(I,3) = 0.0D0
                                                        END IF
                                                END IF 
                                        END IF
!END IF 
                                        BETA(I) = 1.0D0*ACOS(PZ(I))
                                                DUMMY1 = PZ1(I)/(SQRT(PX1(I)**2+PY1(I)**2)*GAMMA(I)**2)
                                        DB(I,1) = PX1(I)*DUMMY1
                                        DB(I,2) = PY1(I)*DUMMY1
                                        DB(I,3) = -1.0D0*SQRT(PX1(I)**2+PY1(I)**2)/GAMMA(I)**2
        
                                        SA(I) = SIN(ALPHA(I))
                                        SB(I) = SIN(BETA(I))
                                        SG(I) = SIN(GAMMA(I))
                                        CA(I) = COS(ALPHA(I))
                                        CB(I) = COS(BETA(I))
                                        CG(I) = COS(GAMMA(I))

                                        QX(I) = -SG(I)*CB(I)*CA(I)-CG(I)*SA(I)
                                        QY(I) = CG(I)*CA(I)-SG(I)*CB(I)*SA(I)
                                        QZ(I) = SG(I)*SB(I)
                                        RX(I) = CG(I)*CB(I)*CA(I)-SG(I)*SA(I)
                                        RY(I) = CG(I)*CB(I)*SA(I)+SG(I)*CA(I)
                                        RZ(I) = -CG(I)*SB(I)


                                        DO J=1,3
                                                DQX(I,J) = -CG(I)*(CB(I)*CA(I))*DG(I,J)-SG(I)*(-SB(I)*CA(I)*DB(I,J)-&
                                                          & CB(I)*SA(I)*DA(I,J))-(-SG(I)*DG(I,J)*SA(I)+CG(I)*CA(I)*DA(I,J))
                                                DQY(I,J) = -SG(I)*DG(I,J)*CA(I)-CG(I)*SA(I)*DA(I,J)-(CG(I)*CB(I)*SA(I)*DG(I,J)+&
                                                          & SG(I)*(-SB(I)*DB(I,J)*SA(I)+CB(I)*CA(I)*DA(I,J)))
                                                DQZ(I,J) = CG(I)*SB(I)*DG(I,J)+SG(I)*CB(I)*DB(I,J)
                  DRX(I,J) = -SG(I)*DG(I,J)*CB(I)*CA(I)+CG(I)*(-SB(I)*CA(I)*DB(I,J)-CB(I)*SA(I)*DA(I,J))-&
                                                          & (CG(I)*SA(I)*DG(I,J)+SG(I)*CA(I)*DA(I,J))
                  DRY(I,J) = -SG(I)*DG(I,J)*CB(I)*SA(I)+CG(I)*(-SB(I)*SA(I)*DB(I,J)+CB(I)*CA(I)*DA(I,J))+&
                                                          &  CG(I)*CA(I)*DG(I,J)-SG(I)*SA(I)*DA(I,J)
                                                DRZ(I,J) = SG(I)*SB(I)*DG(I,J)-CG(I)*CB(I)*DB(I,J)
                                        END DO ! J=1,3
                                END DO  ! I=1,2

! EXTRA LJ SITE AT THE END OF THE REPULSIVE SEMIAXIS C (UNIT VECTOR POINTING TO IT: (RX,RY,RZ))
                                        IF(LJSITE.AND.K==1) THEN
                                             ! FIRST PARTICLE
                                                XLJ(1,1)=X1+PARAMC1*RX(1)
                                                XLJ(1,2)=Y1+PARAMC1*RY(1)
                                                XLJ(1,3)=Z1+PARAMC1*RZ(1)
                                             ! SECOND PARTICLE
                                                XLJ(2,1)=X2+PARAMC1*RX(2)
                                                XLJ(2,2)=Y2+PARAMC1*RY(2)
                                                XLJ(2,3)=Z2+PARAMC1*RZ(2)
                                             DO I=1,3
                                                RLJ2=(XLJ(2,1)-XLJ(1,1))**2+(XLJ(2,2)-XLJ(1,2))**2+(XLJ(2,3)-XLJ(1,3))**2
                                                RLJ=SQRT(RLJ2)
                                             END DO
                                                DRLJ(1)=-1.0D0/RLJ*(XLJ(2,1)-XLJ(1,1))         !DRLJ/DX1
                                                DRLJ(2)=-1.0D0/RLJ*(XLJ(2,2)-XLJ(1,2))         !DRLJ/DY1
                                                DRLJ(3)=-1.0D0/RLJ*(XLJ(2,3)-XLJ(1,3))         !DRLJ/DZ1
                                                DRLJ(4)=-DRLJ(1)                               !DRLJ/DX2
                                                DRLJ(5)=-DRLJ(2)                               !DRLJ/DY2
                                                DRLJ(6)=-DRLJ(3)                               !DRLJ/DZ2

                                                DRLJ(7) = 1.0D0/RLJ*( (XLJ(2,1)-XLJ(1,1))*PARAMC1*(-DRX(1,1)) + &
                                                         &            (XLJ(2,2)-XLJ(1,2))*PARAMC1*(-DRY(1,1)) + &
                                                         &            (XLJ(2,3)-XLJ(1,3))*PARAMC1*(-DRZ(1,1))   ) !DRLJ/DPX1
                                                DRLJ(8) = 1.0D0/RLJ*( (XLJ(2,1)-XLJ(1,1))*PARAMC1*(-DRX(1,2)) + &
                                                         &            (XLJ(2,2)-XLJ(1,2))*PARAMC1*(-DRY(1,2)) + &
                                                         &            (XLJ(2,3)-XLJ(1,3))*PARAMC1*(-DRZ(1,2))   ) !DRLJ/DPY1
                                                DRLJ(9) = 1.0D0/RLJ*( (XLJ(2,1)-XLJ(1,1))*PARAMC1*(-DRX(1,3)) + &
                                                         &            (XLJ(2,2)-XLJ(1,2))*PARAMC1*(-DRY(1,3)) + &
                                                         &            (XLJ(2,3)-XLJ(1,3))*PARAMC1*(-DRZ(1,3))   ) !DRLJ/DPZ1
                                                DRLJ(10)= 1.0D0/RLJ*( (XLJ(2,1)-XLJ(1,1))*PARAMC1*( DRX(2,1)) + &
                                                         &            (XLJ(2,2)-XLJ(1,2))*PARAMC1*( DRY(2,1)) + &
                                                         &            (XLJ(2,3)-XLJ(1,3))*PARAMC1*( DRZ(2,1))   ) !DRLJ/DPX2
                                                DRLJ(11)= 1.0D0/RLJ*( (XLJ(2,1)-XLJ(1,1))*PARAMC1*( DRX(2,2)) + &
                                                         &            (XLJ(2,2)-XLJ(1,2))*PARAMC1*( DRY(2,2)) + &
                                                         &            (XLJ(2,3)-XLJ(1,3))*PARAMC1*( DRZ(2,2))   ) !DRLJ/DPY2
                                                DRLJ(12)= 1.0D0/RLJ*( (XLJ(2,1)-XLJ(1,1))*PARAMC1*( DRX(2,3)) + &
                                                         &            (XLJ(2,2)-XLJ(1,2))*PARAMC1*( DRY(2,3)) + &
                                                         &            (XLJ(2,3)-XLJ(1,3))*PARAMC1*( DRZ(2,3))   ) !DRLJ/DPZ2

                                        END IF

                        A11(1,1)=PX(1)**2*A1(K)+QX(1)**2*B1(K)+RX(1)**2*C1(K)
                        A11(1,2)=PY(1)*PX(1)*A1(K)+QY(1)*QX(1)*B1(K)+RY(1)*RX(1)*C1(K)
                        A11(1,3)=PX(1)*PZ(1)*A1(K)+QX(1)*QZ(1)*B1(K)+RX(1)*RZ(1)*C1(K)
                        A11(2,1)=A11(1,2)
                        A11(2,2)=PY(1)**2*A1(K)+QY(1)**2*B1(K)+RY(1)**2*C1(K)
                        A11(2,3)=PY(1)*PZ(1)*A1(K)+QY(1)*QZ(1)*B1(K)+RY(1)*RZ(1)*C1(K)
                        A11(3,1)=A11(1,3)
                        A11(3,2)=A11(2,3)
                        A11(3,3)=PZ(1)**2*A1(K)+QZ(1)**2*B1(K)+RZ(1)**2*C1(K)

                        B11(1,1)=PX(2)**2*A1(K)+QX(2)**2*B1(K)+RX(2)**2*C1(K)
                        B11(1,2)=PY(2)*PX(2)*A1(K)+QY(2)*QX(2)*B1(K)+RY(2)*RX(2)*C1(K)
                        B11(1,3)=PX(2)*PZ(2)*A1(K)+QX(2)*QZ(2)*B1(K)+RX(2)*RZ(2)*C1(K)
                        B11(2,1)=B11(1,2)
                        B11(2,2)=PY(2)**2*A1(K)+QY(2)**2*B1(K)+RY(2)**2*C1(K)
                        B11(2,3)=PY(2)*PZ(2)*A1(K)+QY(2)*QZ(2)*B1(K)+RY(2)*RZ(2)*C1(K)
                        B11(3,1)=B11(1,3)
                        B11(3,2)=B11(2,3)
                        B11(3,3)=PZ(2)**2*A1(K)+QZ(2)**2*B1(K)+RZ(2)**2*C1(K)

                        DO
                                RANDOM=DPRAND()
                                 DF_LAMBDA(:)=0.0D0
                                       W(:)=0.0D0

                                C_1(1)=(1.0D0-LAMBDA(K))*A11(1,1)+LAMBDA(K)*B11(1,1)
                                C_1(2)=(1.0D0-LAMBDA(K))*A11(1,2)+LAMBDA(K)*B11(1,2)
                                C_1(3)=(1.0D0-LAMBDA(K))*A11(1,3)+LAMBDA(K)*B11(1,3)
                                C_1(4)=(1.0D0-LAMBDA(K))*A11(2,2)+LAMBDA(K)*B11(2,2)
                                        C_1(5)=(1.0D0-LAMBDA(K))*A11(2,3)+LAMBDA(K)*B11(2,3)
                                        C_1(6)=(1.0D0-LAMBDA(K))*A11(3,3)+LAMBDA(K)*B11(3,3)

                                CALL SVERT(C_1,3,W)  !INVERTING THE SYMMETRIC MATRIX C1 (OR G, EQ.2.6 IN ECF PAPER)
     
                                C1INV(1,1)=C_1(1)
                                C1INV(1,2)=C_1(2)
                                C1INV(1,3)=C_1(3)
                                C1INV(2,1)=C1INV(1,2)
                                C1INV(2,2)=C_1(4)
                                C1INV(2,3)=C_1(5)
                                C1INV(3,1)=C1INV(1,3)
                                C1INV(3,2)=C1INV(2,3)
                                C1INV(3,3)=C_1(6)

                                G_INVR=MATMUL(C1INV,R1)
     
                                G_INVRT(1,1)=G_INVR(1)
                                G_INVRT(1,2)=G_INVR(2)
                                G_INVRT(1,3)=G_INVR(3)
     
                                DF_LAMBDA=MATMUL(G_INVRT,MATMUL((1.0D0-LAMBDA(K))**2*A11-LAMBDA(K)**2*B11,G_INVR)) !DERIVATIVE OF THE CONTACT FUNCTION WITH RESPECT TO LAMBDA
                                DUMMY=DF_LAMBDA(1)
                                DUMMY1=LAMBDA(K)

                                IF(ABS(DF_LAMBDA(1))<CONVG) THEN
                                        EXIT
                                ELSE
                                        ITER=ITER+1
!                                        IF THE NEW DERIVATIVE IS LARGER THAN THE OLD ONE, RESET LAMBDA TO ITS PREVIOUS VALUE AND DECREASE STEPSIZE
                                        IF(ABS(DF_LAMBDA(1))>ABS(DF_LAMBDAOLD)) THEN
!                                                WRITE(*,*) 'DF_LAMBDA>DF_LAMBDAOLD', DF_LAMBDA, DF_LAMBDAOLD
                                                LAMBDA(K)=LAMBDAOLD
                                                DF_LAMBDA(1)=DF_LAMBDAOLD
!                                                IF(STEPSIZE>=0.000000001D0)
                                                STEPSIZE=STEPSIZE/3.0D0
                                        END IF

                                        IF(LAMBDA(K)<0.OR.LAMBDA(K)>1) THEN
!                                                WRITE(*,*) "LAMBDA VALUE NOT IN [0,1], RESETTING",J1,J3
                                                STEPSIZE=1.0D0
                                                LAMBDA(K)=0.5
                                                DF_LAMBDA(1)=0.1D0
                                                DUMMY=DF_LAMBDA(1)
                                                DUMMY1=LAMBDA(K)
                                        END IF

                                        IF(ABS(DF_LAMBDA(1))<=0.1D0) THEN
                                                LAMBDA(K)=LAMBDA(K)+DF_LAMBDA(1)*STEPSIZE*RANDOM
                                        ELSE IF(ABS(DF_LAMBDA(1))>1.0D0) THEN
                                                        LAMBDA(K)=LAMBDA(K)+(DF_LAMBDA(1)/ABS(DF_LAMBDA(1)))*RANDOM*STEPSIZE*0.2
                                        ELSE
                                                        LAMBDA(K)=LAMBDA(K)+DF_LAMBDA(1)*STEPSIZE*0.1*RANDOM
                                        END IF

                                        LAMBDAOLD=DUMMY1
                                        DF_LAMBDAOLD=DUMMY
 
                                END IF
!                                IF(ITER.GT.1000) THEN
!                                        WRITE(*,*) LAMBDA(K),DF_LAMBDA(1),ITER,PX,PY,PY/PX,ALPHA
!                                        WRITE(*,*) X1,Y1,Z1
!                                        WRITE(*,*) X2,Y2,Z2
!                                        WRITE(*,*) PX1,PY1,PZ1
!                                        WRITE(*,*) PX2,PY2,PZ2
!                                        WRITE(*,*) SQRT(PX1**2+PY1**2+PZ1**2)
!                                        WRITE(*,*) SQRT(PX2**2+PY2**2+PZ2**2)
!                                END IF
                        END DO

! CALCULATE S(LAMBDA) AT THE CUTOFF, ASSUMING LAMBDA=LAMBDA_C

                        S_LAMBDA=LAMBDA(K)*(1.0D0-LAMBDA(K))*MATMUL(R1T,G_INVR)  ! THIS IS THE ELLIPSOID CONTACT FUNCTION
                        F1(K)=S_LAMBDA(1)
!                        DERIVATIVES OF THE ELLIPSOID CONTACT FUNCTION

!                        WITH RESPECT TO CARTESIAN COORDINATES:
                        DUMMY=2.0D0*LAMBDA(K)*(1.0D0-LAMBDA(K))
                            DO J=1,3
                                DF1(K,J)=DUMMY*G_INVR(J)*(-1.0D0)  !DF/DX1,DY1,DZ1
                                DF1(K,J+3)=-DF1(K,J)           !DF/DX2,DY2,DZ2
                            END DO

!                        WITH RESPECT TO ORIENTATION:
                        DO J=1,2
                                DPX(J,1) = (PY1(J)**2+PZ1(J)**2)/GAMMA(J)**3                      ! DPX/DPX1
                                DPX(J,2) = -1.0D0*(PX1(J)*PY1(J))/GAMMA(J)**3                     ! DPX/DPY1
                                DPX(J,3) = -1.0D0*(PX1(J)*PZ1(J))/GAMMA(J)**3                     ! DPX/DPZ1
                                DPY(J,1) = DPX(J,2)                                        ! DPY/DPX1 = DPX/DPY1
                                DPY(J,2) = (PX1(J)**2+PZ1(J)**2)/GAMMA(J)**3                      ! DPY/DPY1
                                DPY(J,3) = -1.0D0*(PY1(J)*PZ1(J))/GAMMA(J)**3                     ! DPY/DPZ1
                                DPZ(J,1) = DPX(J,3)                                        ! DPZ/DPX1 = DPX/DPZ1
                                DPZ(J,2) = DPY(J,3)                                        ! DPZ/DPY1 = DPY/DPZ1
                                DPZ(J,3) = (PX1(J)**2+PY1(J)**2)/GAMMA(J)**3                      ! DPZ/DPZ1
                        END DO

                        DO I=1,3
                                DAA(I,K,1,1)=2.0D0*PX(1)*A1(K)*DPX(1,I)+2.0D0*QX(1)*B1(K)*DQX(1,I)+&
                                            &2.0D0*RX(1)*C1(K)*DRX(1,I)
                                DAA(I,K,1,2)=A1(K)*(DPY(1,I)*PX(1)+PY(1)*DPX(1,I))+B1(K)*(DQY(1,I)*QX(1)+&
                                            &QY(1)*DQX(1,I))+C1(K)*(DRY(1,I)*RX(1)+RY(1)*DRX(1,I))
                                DAA(I,K,1,3)=A1(K)*(DPZ(1,I)*PX(1)+PZ(1)*DPX(1,I))+B1(K)*(DQZ(1,I)*QX(1)+&
                                            &QZ(1)*DQX(1,I))+C1(K)*(DRZ(1,I)*RX(1)+RZ(1)*DRX(1,I))
                                DAA(I,K,2,1)=DAA(I,K,1,2)
                                DAA(I,K,2,2)=2.0D0*PY(1)*A1(K)*DPY(1,I)+2.0D0*QY(1)*B1(K)*DQY(1,I)+&
                                            &2.0D0*RY(1)*C1(K)*DRY(1,I)
                                DAA(I,K,2,3)=A1(K)*(DPY(1,I)*PZ(1)+PY(1)*DPZ(1,I))+B1(K)*(DQY(1,I)*QZ(1)+&
                                            &QY(1)*DQZ(1,I))+C1(K)*(DRY(1,I)*RZ(1)+RY(1)*DRZ(1,I))
                                DAA(I,K,3,1)=DAA(I,K,1,3)
                                DAA(I,K,3,2)=DAA(I,K,2,3)
                                DAA(I,K,3,3)=2.0D0*PZ(1)*A1(K)*DPZ(1,I)+2.0D0*QZ(1)*B1(K)*DQZ(1,I)+&
                                            &2.0D0*RZ(1)*C1(K)*DRZ(1,I)
                                DBB(I,K,1,1)=2.0D0*PX(2)*A1(K)*DPX(2,I)+2.0D0*QX(2)*B1(K)*DQX(2,I)+&
                                            &2.0D0*RX(2)*C1(K)*DRX(2,I)
                                DBB(I,K,1,2)=A1(K)*(DPY(2,I)*PX(2)+PY(2)*DPX(2,I))+B1(K)*(DQY(2,I)*QX(2)+&
                                            &QY(2)*DQX(2,I))+C1(K)*(DRY(2,I)*RX(2)+RY(2)*DRX(2,I))
                                DBB(I,K,1,3)=A1(K)*(DPZ(2,I)*PX(2)+PZ(2)*DPX(2,I))+B1(K)*(DQZ(2,I)*QX(2)+&
                                            &QZ(2)*DQX(2,I))+C1(K)*(DRZ(2,I)*RX(2)+RZ(2)*DRX(2,I))
                                DBB(I,K,2,1)=DBB(I,K,1,2)
                                DBB(I,K,2,2)=2.0D0*PY(2)*A1(K)*DPY(2,I)+2.0D0*QY(2)*B1(K)*DQY(2,I)+&
                                            &2.0D0*RY(2)*C1(K)*DRY(2,I)
                                DBB(I,K,2,3)=A1(K)*(DPY(2,I)*PZ(2)+PY(2)*DPZ(2,I))+B1(K)*(DQY(2,I)*QZ(2)+&
                                            &QY(2)*DQZ(2,I))+C1(K)*(DRY(2,I)*RZ(2)+RY(2)*DRZ(2,I))
                                DBB(I,K,3,1)=DBB(I,K,1,3)
                                DBB(I,K,3,2)=DBB(I,K,2,3)
                                DBB(I,K,3,3)=2.0D0*PZ(2)*A1(K)*DPZ(2,I)+2.0D0*QZ(2)*B1(K)*DQZ(2,I)+&
                                            &2.0D0*RZ(2)*C1(K)*DRZ(2,I)
                        END DO

                        DO I=1,3
                                DGU1(I,:,:)=-(1.0D0-LAMBDA(K))*MATMUL(MATMUL(C1INV, DAA(I,K,:,:)),C1INV)
                                DGU2(I,:,:)=-LAMBDA(K)*MATMUL(MATMUL(C1INV, DBB(I,K,:,:)),C1INV)
                        END DO     
            
                        DO J=1,3
                                DXU1(J,:)=MATMUL(DGU1(J,:,:),R1)
                                DFU1(J,:)=MATMUL(R1T,DXU1(J,:))
                                DXU2(J,:)=MATMUL(DGU2(J,:,:),R1)
                                DFU2(J,:)=MATMUL(R1T,DXU2(J,:))
                        END DO
                    
                            DUMMY=LAMBDA(K)*(1.0D0-LAMBDA(K))
                    
!                           DF1(K,7 )=DUMMY*DFU1(1,1)   !DF/DALPHA1
!                           DF1(K,8 )=DUMMY*DFU1(2,1)      !DF/DBETA1
!                           DF1(K,9 )=DUMMY*DFU2(1,1)       !DF/DALPHA2
!                           DF1(K,10)=DUMMY*DFU2(2,1)   !DF/DBETA2
                    
                        DF1(K,7 )=DUMMY*DFU1(1,1)      !DF/D(PX1)
                        DF1(K,8 )=DUMMY*DFU1(2,1)      !DF/D(PY1)
                        DF1(K,9 )=DUMMY*DFU1(3,1)      !DF/D(PZ1)

                        DF1(K,10)=DUMMY*DFU2(1,1)       !DF/D(PX2)
                        DF1(K,11)=DUMMY*DFU2(2,1)       !DF/D(PY2)
                        DF1(K,12)=DUMMY*DFU2(3,1)       !DF/D(PZ2)
                    
!                        END OF DERIVATIVES OF THE ELLIPSOID CONTACT FUNCTION

         
                            LJ1(K)=SIGMA0(K)/(R-R/SQRT(F1(K))+SIGMA0(K))
                            LJ2(K)=LJ1(K)**2
                            LJ3(K)=LJ2(K)*LJ1(K)
                            LJ4(K)=LJ2(K)**2
                            LJ5(K)=LJ4(K)*LJ1(K)
                            LJ6(K)=LJ4(K)*LJ2(K)
                            LJ7(K)=LJ6(K)*LJ1(K)
                            LJ11(K)=LJ5(K)*LJ6(K)
                            LJ12(K)=LJ6(K)**2

! CORRECTION TERMS TO THE POTENTIAL IF WE REQUIRE A CUTOFF AT RC. 
                            IF (PARAMONOVCUTOFF) THEN
! R/SQRT(F1(K)) = F(A,B), A PARAMETER ONLY DEPENDENT ON ORIENTATION!

                                       CLJ1(K)=SIGMA0(K)/(PCUTOFF-R/SQRT(F1(K))+SIGMA0(K))
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
                            END IF
 
!                        LLJ1=SIGMA1/(R-R/SQRT(F1(1))+SIGMA1)
!                        LLJ2=LLJ1**2
!                        LLJ3=LLJ2*LLJ1
!                        LLJ4=LLJ2**2
!                        LLJ5=LLJ4*LLJ1
!                        LLJ6=LLJ4*LLJ2
!                        LLJ7=LLJ6*LLJ1
!                        LLJ11=LLJ5*LLJ6
!                        LLJ12=LLJ6**2
              ! INTERACTION BETWEEN THE EXTRA LJ SITES:
                        LLJ1=SIGMA1/RLJ
                        LLJ2=LLJ1**2
                        LLJ3=LLJ2*LLJ1
                        LLJ4=LLJ2**2
                        LLJ5=LLJ4*LLJ1
                        LLJ6=LLJ4*LLJ2
                        LLJ7=LLJ6*LLJ1
                        LLJ11=LLJ5*LLJ6
                        LLJ12=LLJ6**2
 
                            DUMMY=1.0D0/RLJ
                            DUMMY=DUMMY**2
                            DO J=1,12
                                DLLJ1(J) =-SIGMA1*DUMMY*DRLJ(J)
                            END DO


!                           WRITE(*,*) "DF(LAMBDA(K))=",LAMBDA(K),DF_LAMBDA(K),ITER

                            DUMMY=LJ1(K)/SIGMA0(K)
                            DUMMY=DUMMY**2
                            DUMMY1=1.0D0/SQRT(F1(K))
                            DUMMY2=DUMMY1**3
 
                            DO J=1,6
                                DLJ1(K,J) =-1.0D0*SIGMA0(K)*DUMMY*(DR(J)-DR(J)*DUMMY1+R*0.5D0*DUMMY2*DF1(K,J))
                            END DO
        
                            DO J=7,12
                                DLJ1(K,J) =-1.0D0*SIGMA0(K)*DUMMY*(0.5D0*R*DUMMY2*DF1(K,J))
                            END DO
! INNER HARD CORE STUFF
!                            DO J=1,6
!                                DLLJ1(J) =-1.0D0*SIGMA1*DUMMY*(DR(J)-DR(J)*DUMMY1+R*0.5D0*DUMMY2*DF1(1,J))
!                        END DO
        
!                        DO J=7,12
!                                DLLJ1(J) =-1.0D0*SIGMA1*DUMMY*(0.5D0*R*DUMMY2*DF1(1,J))
!                        END DO


                        IF (PARAMONOVCUTOFF) THEN 
                                DUMMY=CLJ1(K)/SIGMA0(K)
                                DUMMY=DUMMY**2
                                DUMMY1=1.0D0/SQRT(F1(K))
                                DUMMY2=DUMMY1**3
 
                                DO J=1,6
                                        DCLJ1(K,J) =-1.0D0*SIGMA0(K)*DUMMY*(-1.0D0*DUMMY1*DR(J)+0.5D0*R*DUMMY2*DF1(K,J))
                                END DO
                                DO J=7,12
                                        DCLJ1(K,J) =-1.0D0*SIGMA0(K)*DUMMY*(0.5D0*R*DUMMY2*DF1(K,J)) !DERIVATIVES WRT TO ORIENTATION
                                END DO
                        END IF
!                        WRITE(*,*) LAMBDA(K),K,ITER
                  END DO        


!                  WRITE(*,*) 'PARAMONOV.F90> ECF = ', F1(1), F1(2)    
                END IF                         ! IF (NOT(PARAMONOVCUTOFF) OR (PARAMONOVCUTOFF AND R<PCUTOFF)) THEN....

        IF (PARAMONOVCUTOFF) THEN
                 IF (R>=PCUTOFF) THEN
                         VDUMMY = 0.0D0
                 ELSE

                  ! WORK OUT THE SPLINE TERMS FOR SMOOTH CUTOFF OF EXTRA LJ SITES
                   IF(RLJ2<RON2) THEN

                      TERM2=1.0D0
                      TERM3=0.0D0
                   ELSE IF(RLJ2>RON2) THEN
                          TERM2=(CUT-RLJ**2)*(CUT-RLJ**2)*(CUT+2.0D0*RLJ**2-3.0D0*RON2)*RANGE2INV3
                          TERM3=RLJ*12.0D0*(CUT-RLJ**2)*(RON2-RLJ**2)*RANGE2INV3 ! D(TERM2)/DR
                   END IF

                   VDUMMY = 0.0D0

!         ADD CORRECTIONS TO NORMAL LJ PAIR POTENTIAL TO ENSURE THAT POTENTIAL AND 
!         DERIVATIVES ARE ZERO AT RC, AND VANISH SMOOTHLY WITH NO DISCONTINUITIES. 
!                VDUMMY=EPSILON1*(LLJ12-0.0D0*LLJ6) !INNER HARD CORE

                VDUMMY=EPSILON1*(LLJ12-LLJ6)*TERM2
!REPULSIVE POTENTIAL AND PERIODIC CUTOFF CORRECTIONS
                VDUMMY=VDUMMY+4.0D0*EPSILON0*(LJ12(1)+(6.0D0*CLJ12(1))*CLJ2(1)/LJ2(1)-7.0D0*CLJ12(1))
!ATTRACTIVE POTENTIAL AND PERIODIC CUTOFF CORRECTIONS
                VDUMMY=VDUMMY+4.0D0*EPSILON0*(-LJ6(2)+(-3.0D0*CLJ6(2))*CLJ2(2)/LJ2(2) +4.0D0*CLJ6(2))

                        END IF
        ELSE
!                VDUMMY=(F12 - F6)
!                VDUMMY=EPSILON1*(LLJ12-LLJ6)*TERM2
!                WRITE(*,*) 'VDUMMY FOR LJ SITES',VDUMMY,RLJ
                VDUMMY=VDUMMY+4.0D0*EPSILON0*(LJ12(1)-LJ6(2))    
        END IF
                      DO J=1,12
                        IF (PARAMONOVCUTOFF) THEN
                                IF (R>=PCUTOFF) THEN
                                        DVDUMMY(J)=0.0D0
                                ELSE
                                        DVDUMMY(J)=0.0D0

                                        DVDUMMY(J)=EPSILON1*(12.0D0*LLJ11*DLLJ1(J)-6.0D0*LLJ5*DLLJ1(J))*TERM2 + &
                                                        & EPSILON1*(LLJ12-LLJ6)*TERM3*DRLJ(J)! EXTRA LJ SITE DERIVATIVES
                                        DVDUMMY(J)=DVDUMMY(J)+4.0D0*EPSILON0*(&
                                                &12.0D0*DLJ1(1,J)*(LJ11(1)-CLJ14(1)/LJ3(1))&
                                                &+84.0D0*DCLJ1(1,J)*(CLJ13(1)/LJ2(1)-CLJ11(1))) !REPULSIVE DERIVATIVE
                                        DVDUMMY(J)=DVDUMMY(J)+4.0D0*EPSILON0*(&
                                                &-6.0D0*DLJ1(2,J)*(LJ5(2)-CLJ8(2)/LJ3(2))&
                                                &-24.0D0*DCLJ1(2,J)*(CLJ7(2)/LJ2(2)-CLJ5(2))) !ATTRACTIVE DERIVATIVE


!4.0D0*EPSILON0*(-6.0D0*DLJ1(2,J)*LJ5(2)+24.0D0*DCLJ1(2,J)*CLJ5(2)-24.0D0*DCLJ1(2,J)*CLJ7(2)/LJ2(2)+6.0D0*DLJ1(2,J)*CLJ8(2)/LJ3(2)) 
!-CLJ8(2)/LJ3(2))-24.0D0*DCLJ1(2,J)*(CLJ7(2)/LJ2(2)))

                                END IF
                        ELSE
!                                 DVDUMMY(J)=EPSILON1*(12.0D0*LLJ11*DLLJ1(J)-6.0D0*LLJ5*DLLJ1(J))&
!                                              &+4.0D0*EPSILON0*(12.0D0*LJ11(1)*DLJ1(1,J)-6.0D0*LJ5(2)*DLJ1(2,J))
                                  DVDUMMY(J)=4.0D0*EPSILON0*(12.0D0*LJ11(1)*DLJ1(1,J)-6.0D0*LJ5(2)*DLJ1(2,J))

                        END IF
!                WRITE(*,*) J, DVDUMMY(J)
                END DO !J=1,12
!       WRITE(*,*) 'I,J,DVDUMMY', J1,J3,ALPHA(2),DA(2,:)
!       WRITE(*,*) PZ1(2),DVDUMMY(1:12)
!                END DO
!                WRITE(*,*) DCLJ1(1,:),DCLJ1(2,:)

!                        Z1 = Z1+0.0000001D0
!                        FCHECK(L)=VDUMMY
!                        IF(L==2) DCHECK=DVDUMMY(3)
!                        IF(L==3) THEN
!                         WRITE(*,*) (FCHECK(3)-FCHECK(1))/0.0000002D0, DCHECK
!                         STOP
!                        END IF
!             END DO    ! FOR DERIVATIVE CHECK


                EGB=EGB+VDUMMY
    
                VT(J1) = VT(J1) + VDUMMY
                VT(J3) = VT(J3) + VDUMMY        ! PAIR POTENTIALS

                V(J2-2)=V(J2-2)+DVDUMMY(1)
                V(K2-2)=V(K2-2)+DVDUMMY(4)
                V(J2-1)=V(J2-1)+DVDUMMY(2)
                V(K2-1)=V(K2-1)+DVDUMMY(5)
                V(J2  )=V(J2  )+DVDUMMY(3)        
                V(K2  )=V(K2  )+DVDUMMY(6)       
                V(3*REALNATOMS+J2-2)=V(3*REALNATOMS+J2-2)+DVDUMMY(7)        
                V(3*REALNATOMS+K2-2)=V(3*REALNATOMS+K2-2)+DVDUMMY(10)        
                V(3*REALNATOMS+J2-1)=V(3*REALNATOMS+J2-1)+DVDUMMY(8)        
                V(3*REALNATOMS+K2-1)=V(3*REALNATOMS+K2-1)+DVDUMMY(11)       
                V(3*REALNATOMS+J2  )=V(3*REALNATOMS+J2  )+DVDUMMY(9)        
                V(3*REALNATOMS+K2  )=V(3*REALNATOMS+K2  )+DVDUMMY(12)       

!                WRITE(*,*) 'DUMPING DERIVATIVES:',J1,J3
!                WRITE(*,*) DVDUMMY(:)
!                WRITE(*,*) X(:)
!                =SQRT((X1-X2)**2,(Y1-Y2)**2,(Z1-Z2)**2)

! SOME DEBUG PRINTING

!WRITE(*,*) 'X1,Y1,Z1'
!WRITE(*,*) X1,Y1,Z1
!WRITE(*,*) 'X2,Y2,Z2'
!WRITE(*,*) X2,Y2,Z2
!WRITE(*,*) 'PX1,PY1,PZ1'
!WRITE(*,*) PX1(1),PY1(1),PZ1(1)
!WRITE(*,*) 'PX2,PY2,PZ2'
!WRITE(*,*) PX1(2),PY1(2),PZ1(2)

!WRITE(*,*) 'ENERGY:' , VDUMMY
!WRITE(*,*) 'LAMBDA, REPULSIVE CONTACT FUNCTION', LAMBDA(1), F1(1) 
!WRITE(*,*) 'LAMBDA, ATTRACTIVE CONTACT FUNCTION', LAMBDA(2), F1(2)
!STOP
            END DO      !DO J3=J1+1,REALNATOMS
           END DO       !DO J1=1,REALNATOMS
        
!           IF(STEST) CALL PARAMONOVSECDER(X,.TRUE.)

!      WRITE(*,*) "EXITING ENERGY CYCLE"
999    CONTINUE      
       END SUBROUTINE OLDPARAMONOV


SUBROUTINE PARAMONOVNUMFIRSTDER(OLDX,STEST)
USE COMMONS, ONLY : NATOMS
!USE KEY
USE MODHESS
IMPLICIT NONE

DOUBLE PRECISION  :: V(3*NATOMS),ENERGY,X(3*NATOMS),NUMGRAD(3*NATOMS),TRUEGRAD(3*NATOMS),&
        &            OLDX(3*NATOMS),ETEMP(2,3*NATOMS),KSI
INTEGER    :: I,J
LOGICAL    :: GTEST,STEST


KSI=0.0001D0
X(:)=OLDX(:)

!         CALL OLDPARAMONOV(X,TRUEGRAD,ENERGY,GTEST,.FALSE.)
         CALL PYGPERIODIC(X,TRUEGRAD,ENERGY,.TRUE.)


!WRITE(*,*) 'SF344> IN PARAMONOVSECDER'
DO I=1,3*NATOMS

         X(I)=X(I)-KSI

!         CALL OLDPARAMONOV(X,V,ENERGY,GTEST,.FALSE.)
         CALL PYGPERIODIC(X,V,ENERGY,.TRUE.)

         ETEMP(1,I)=ENERGY

         X(I)=X(I)+2.0D0*KSI

!         CALL OLDPARAMONOV(X,V,ENERGY,GTEST,.FALSE.)
         CALL PYGPERIODIC(X,V,ENERGY,.TRUE.)
    
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





SUBROUTINE ECF(OVERLAPTEST,ECFVALUE,X1,X2,Y1,Y2,Z1,Z2,PX11,PX12,PY11,PY12,PZ11,PZ12,A,B,C)


!CALCULATES THE VALUE OF THE ELLIPSOID CONTACT FUNCTION

!USE COMMONS
IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN)       :: X1,X2,Y1,Y2,Z1,Z2,PX11,PX12,PY11,PY12,PZ11,PZ12,A,B,C
DOUBLE PRECISION, INTENT(OUT)      :: ECFVALUE
DOUBLE PRECISION                   :: A1(2),B1(2),C1(2),LAMBDAOLD,DF_LAMBDAOLD,&
                &       DUMMY,DUMMY1,DUMMY2,S_DUMMY1,C_DUMMY1,S_DUMMY2,C_DUMMY2,PX1(2),PY1(2),PZ1(2),&
                &       ALPHA(2),BETA(2),GAMMA(2),QX(2),QY(2),QZ(2),RX(2),RY(2),RZ(2),PX(2),PY(2),PZ(2),&
                &       SA(2),CA(2),SB(2),CB(2),SG(2),CG(2)
DOUBLE PRECISION  :: S_LAMBDA(1),F1(1),LAMBDA(2), A11(3,3),B11(3,3),C_1(6),C1INV(3,3),W(3),R1(3),R1T(1,3),&
                &       DF_LAMBDA(1),G_INVR(3),G_INVRT(1,3),DPRAND, RANDOM, STEPSIZE, CONVG, R, PI
INTEGER    :: I,J,K,ITER
!LOGICAL, INTENT(OUT)    :: OVERLAPTEST
LOGICAL   :: OVERLAPTEST

A1(1) = A**2
B1(1) = B**2
C1(1) = C**2

PX1(1) = PX11
PX1(2) = PX12
PY1(1) = PY11
PY1(2) = PY12
PZ1(1) = PZ11
PZ1(2) = PZ12

PI = ATAN(1.0)*4.0D0

K=1
        R=SQRT((X2-X1)**2+(Y2-Y1)**2+(Z2-Z1)**2)
                   
 R1(1)=X2-X1
 R1(2)=Y2-Y1
 R1(3)=Z2-Z1
 
 R1T(1,1)=R1(1)
 R1T(1,2)=R1(2)
 R1T(1,3)=R1(3)

           STEPSIZE=1.0D0
            LAMBDA(1)=0.5D0
            ITER=0
            CONVG=1.0D-6
            DF_LAMBDAOLD=1.0D6
 
   DO I=1,2
     GAMMA(I) = SQRT(PX1(I)**2+PY1(I)**2+PZ1(I)**2)

     PX(I) = PX1(I)/GAMMA(I)
     PY(I) = PY1(I)/GAMMA(I)
     PZ(I) = PZ1(I)/GAMMA(I)

     IF(PX(I).EQ.0.0D0) THEN
       IF(PY(I)>=0.0D0) THEN
        ALPHA(I) = PI/2           ! EULER ANGLE ALPHA 
       ELSE
        ALPHA(I) = -1.0D0*PI/2
       END IF
     ELSE
        IF(PY(I)>=0.0D0) THEN
          IF(PX(I)>0.0D0) THEN
            ALPHA(I) = 1.0D0*ATAN(PY(I)/PX(I))   ! FIRST QUADRANT
          ELSE    ! PX<0
            ALPHA(I) = 1.0D0*ATAN(PY(I)/PX(I))+PI       ! SHOULD BE IN THE SECOND QUADRANT
          END IF
        ELSE IF(PY(I)<0.0D0) THEN
          IF(PX(I)>0.0D0) THEN
            ALPHA(I) = 1.0D0*ATAN(PY(I)/PX(I))    ! FOURTH QUADRANT
          ELSE    ! PX<0
            ALPHA(I) = 1.0D0*ATAN(PY(I)/PX(I))-PI            ! THIRD QUADRANT
          END IF
        END IF 
     END IF
 
        BETA(I) = 1.0D0*ACOS(PZ(I))

        DUMMY1 = PZ1(I)/(SQRT(PX1(I)**2+PY1(I)**2)*GAMMA(I)**2)
        
     SA(I) = SIN(ALPHA(I))
     SB(I) = SIN(BETA(I))
     SG(I) = SIN(GAMMA(I))
     CA(I) = COS(ALPHA(I))
     CB(I) = COS(BETA(I))
     CG(I) = COS(GAMMA(I))

     QX(I) = -SG(I)*CB(I)*CA(I)-CG(I)*SA(I)
     QY(I) = CG(I)*CA(I)-SG(I)*CB(I)*SA(I)
     QZ(I) = SG(I)*SB(I)
     RX(I) = CG(I)*CB(I)*CA(I)-SG(I)*SA(I)
     RY(I) = CG(I)*CB(I)*SA(I)+SG(I)*CA(I)
     RZ(I) = -CG(I)*SB(I)
   END DO  ! I=1,2

      A11(1,1)=PX(1)**2*A1(K)+QX(1)**2*B1(K)+RX(1)**2*C1(K)
      A11(1,2)=PY(1)*PX(1)*A1(K)+QY(1)*QX(1)*B1(K)+RY(1)*RX(1)*C1(K)
      A11(1,3)=PX(1)*PZ(1)*A1(K)+QX(1)*QZ(1)*B1(K)+RX(1)*RZ(1)*C1(K)
      A11(2,1)=A11(1,2)
      A11(2,2)=PY(1)**2*A1(K)+QY(1)**2*B1(K)+RY(1)**2*C1(K)
      A11(2,3)=PY(1)*PZ(1)*A1(K)+QY(1)*QZ(1)*B1(K)+RY(1)*RZ(1)*C1(K)
      A11(3,1)=A11(1,3)
      A11(3,2)=A11(2,3)
      A11(3,3)=PZ(1)**2*A1(K)+QZ(1)**2*B1(K)+RZ(1)**2*C1(K)

      B11(1,1)=PX(2)**2*A1(K)+QX(2)**2*B1(K)+RX(2)**2*C1(K)
      B11(1,2)=PY(2)*PX(2)*A1(K)+QY(2)*QX(2)*B1(K)+RY(2)*RX(2)*C1(K)
      B11(1,3)=PX(2)*PZ(2)*A1(K)+QX(2)*QZ(2)*B1(K)+RX(2)*RZ(2)*C1(K)
      B11(2,1)=B11(1,2)
      B11(2,2)=PY(2)**2*A1(K)+QY(2)**2*B1(K)+RY(2)**2*C1(K)
      B11(2,3)=PY(2)*PZ(2)*A1(K)+QY(2)*QZ(2)*B1(K)+RY(2)*RZ(2)*C1(K)
      B11(3,1)=B11(1,3)
      B11(3,2)=B11(2,3)
      B11(3,3)=PZ(2)**2*A1(K)+QZ(2)**2*B1(K)+RZ(2)**2*C1(K)


           DO
                  RANDOM=DPRAND()
                   DF_LAMBDA(:)=0.0D0
                   W(:)=0.0D0

                  C_1(1)=(1.0D0-LAMBDA(K))*A11(1,1)+LAMBDA(K)*B11(1,1)
                  C_1(2)=(1.0D0-LAMBDA(K))*A11(1,2)+LAMBDA(K)*B11(1,2)
                  C_1(3)=(1.0D0-LAMBDA(K))*A11(1,3)+LAMBDA(K)*B11(1,3)
                  C_1(4)=(1.0D0-LAMBDA(K))*A11(2,2)+LAMBDA(K)*B11(2,2)
                  C_1(5)=(1.0D0-LAMBDA(K))*A11(2,3)+LAMBDA(K)*B11(2,3)
                  C_1(6)=(1.0D0-LAMBDA(K))*A11(3,3)+LAMBDA(K)*B11(3,3)

 
 
        ! WRITE(*,*) 'CALLING SVERT, MATRIX TO BE INVERTED:'
        ! WRITE(*,*) C_1(:)
        ! WRITE(*,*) ITER
        ! WRITE(*,*) W(:)
          
             CALL SVERT(C_1,3,W)  !INVERTING THE SYMMETRIC MATRIX C1 (OR G, EQ.2.6 IN ECF PAPER)
     

      
             C1INV(1,1)=C_1(1)
             C1INV(1,2)=C_1(2)
             C1INV(1,3)=C_1(3)
             C1INV(2,1)=C1INV(1,2)
             C1INV(2,2)=C_1(4)
             C1INV(2,3)=C_1(5)
             C1INV(3,1)=C1INV(1,3)
             C1INV(3,2)=C1INV(2,3)
             C1INV(3,3)=C_1(6)


        ! WRITE(*,*) 'INVERTED MATRIX:'
        ! WRITE(*,*) C1INV(:,:)
     
            G_INVR=MATMUL(C1INV,R1)
     
             G_INVRT(1,1)=G_INVR(1)
             G_INVRT(1,2)=G_INVR(2)
             G_INVRT(1,3)=G_INVR(3)
     
    
      !DERIVATIVE OF THE CONTACT FUNCTION WITH RESPECT TO LAMBDA
       DF_LAMBDA=MATMUL(G_INVRT,MATMUL((1.0D0-LAMBDA(K))**2*A11-LAMBDA(K)**2*B11,G_INVR))   
          DUMMY=DF_LAMBDA(1)
          DUMMY1=LAMBDA(K)

             IF(ABS(DF_LAMBDA(1))<CONVG) THEN
                EXIT
             ELSE
                ITER=ITER+1

        ! IF THE NEW DERIVATIVE IS LARGER THAN THE OLD ONE, RESET LAMBDA TO ITS PREVIOUS VALUE AND DECREASE STEPSIZE
                IF(ABS(DF_LAMBDA(1))>ABS(DF_LAMBDAOLD)) THEN
!                WRITE(*,*) 'DF_LAMBDA>DF_LAMBDAOLD', DF_LAMBDA, DF_LAMBDAOLD
                    LAMBDA(K)=LAMBDAOLD
                    DF_LAMBDA(1)=DF_LAMBDAOLD
                    !IF(STEPSIZE>=0.000000001D0)
                    STEPSIZE=STEPSIZE/3.0D0
                END IF

                IF(LAMBDA(K)<0.OR.LAMBDA(K)>1) THEN
                   WRITE(*,*) "LAMBDA VALUE NOT IN [0,1], RESETTING"
                    STEPSIZE=1.0D0
                    LAMBDA(K)=0.5
                    DF_LAMBDA(1)=0.1D0
                    DUMMY=DF_LAMBDA(1)
                    DUMMY1=LAMBDA(K)
                END IF

                IF(ABS(DF_LAMBDA(1))<=0.1D0) THEN
                    LAMBDA(K)=LAMBDA(K)+DF_LAMBDA(1)*STEPSIZE*RANDOM
                ELSE IF(ABS(DF_LAMBDA(1))>1.0D0) THEN
                    LAMBDA(K)=LAMBDA(K)+(DF_LAMBDA(1)/ABS(DF_LAMBDA(1)))*RANDOM*STEPSIZE*0.2
                ELSE
                    LAMBDA(K)=LAMBDA(K)+DF_LAMBDA(1)*STEPSIZE*0.1*RANDOM
                END IF
 
  LAMBDAOLD=DUMMY1
  DF_LAMBDAOLD=DUMMY

             END IF
!           IF(ITER.GT.1000) THEN
!
!                WRITE(*,*) LAMBDA(K),DF_LAMBDA(1),ITER,PX,PY,PY/PX,ALPHA
!                WRITE(*,*) X1,Y1,Z1
!                WRITE(*,*) X2,Y2,Z2
!                WRITE(*,*) PX1,PY1,PZ1
!                WRITE(*,*) PX2,PY2,PZ2
!                WRITE(*,*) SQRT(PX1**2+PY1**2+PZ1**2)
!                WRITE(*,*) SQRT(PX2**2+PY2**2+PZ2**2)
!           END IF
           END DO
   
           S_LAMBDA=LAMBDA(K)*(1.0D0-LAMBDA(K))*MATMUL(R1T,G_INVR)  ! THIS IS THE ELLIPSOID CONTACT FUNCTION
           F1(K)=S_LAMBDA(1)
ECFVALUE=S_LAMBDA(1)

IF(ECFVALUE<1.0D0) THEN
 OVERLAPTEST=.TRUE.
ELSE
 OVERLAPTEST=.FALSE.
END IF

!WRITE(*,*) K
!WRITE(*,*) OVERLAPTEST

END SUBROUTINE ECF


SUBROUTINE OLDECF2(OVERLAPTEST,ECFVALUE,X1,X2,Y1,Y2,Z1,Z2,PX1,PX2,PY1,PY2,PZ1,PZ2,A,B,C)


!CALCULATES THE VALUE OF THE ELLIPSOID CONTACT FUNCTION

!USE COMMONS
IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN)       :: X1,X2,Y1,Y2,Z1,Z2,PX1,PX2,PY1,PY2,PZ1,PZ2,A,B,C
DOUBLE PRECISION, INTENT(OUT)      :: ECFVALUE
DOUBLE PRECISION                   :: A1(2),B1(2),C1(2),LAMBDAOLD,DF_LAMBDAOLD,&
                &       DUMMY,DUMMY1,DUMMY2,S_DUMMY1,C_DUMMY1,S_DUMMY2,C_DUMMY2
DOUBLE PRECISION  :: S_LAMBDA(1),F1(1),LAMBDA(2), A11(3,3),B11(3,3),C_1(6),C1INV(3,3),W(3),R1(3),R1T(1,3),&
                &       DF_LAMBDA(1),G_INVR(3),G_INVRT(1,3),DPRAND, RANDOM, STEPSIZE, CONVG, R, PI
INTEGER    :: I,J,K,ITER
LOGICAL, INTENT(OUT)    :: OVERLAPTEST

A1(1) = A**2
B1(1) = B**2
C1(1) = C**2

PI = ATAN(1.0)*4.0D0

K=1
        R=SQRT((X2-X1)**2+(Y2-Y1)**2+(Z2-Z1)**2)
                   
 R1(1)=X2-X1
 R1(2)=Y2-Y1
 R1(3)=Z2-Z1
 
 R1T(1,1)=R1(1)
 R1T(1,2)=R1(2)
 R1T(1,3)=R1(3)

!       DR(1)=-1.0D0/R*(X2-X1)
!       DR(2)=-1.0D0/R*(Y2-Y1)
!       DR(3)=-1.0D0/R*(Z2-Z1)
!       DR(4)=-DR(1)
!       DR(5)=-DR(2)
!       DR(6)=-DR(3)
                          
!      PX1=C_ALPHA1*C_BETA1
!      PY1=C_ALPHA1*S_BETA1 
!      PZ1=S_ALPHA1        
!      PX2=C_ALPHA2*C_BETA2
!      PY2=C_ALPHA2*S_BETA2
!      PZ2=S_ALPHA2

           STEPSIZE=1.0D0
            LAMBDA(K)=0.5D0
            ITER=0
            CONVG=1.0D-2
            DF_LAMBDAOLD=1.0D6
 

     IF(PX1.EQ.0.0D0) THEN
        DUMMY1 = PI/2
     ELSE
        DUMMY1 = ATAN(PY1/PX1)
     END IF

     IF(PX2.EQ.0.0D0) THEN
        DUMMY2 = PI/2
     ELSE
        DUMMY2 = ATAN(PY2/PX2)
     END IF

        S_DUMMY1 = SIN(DUMMY1)
        C_DUMMY1 = COS(DUMMY1)
        S_DUMMY2 = SIN(DUMMY2)
        C_DUMMY2 = COS(DUMMY2)

   
      DUMMY = SQRT(PX1**2+PY1**2+PZ1**2)

      A11(1,1)=(PX1/DUMMY)**2*A1(K)+S_DUMMY1**2*B1(K)+C_DUMMY1**2*(PZ1/DUMMY)**2*C1(K)
      A11(1,2)=PY1*PX1*A1(K)/(DUMMY**2)-C_DUMMY1*S_DUMMY1*B1(K)+PZ1**2*C_DUMMY1*S_DUMMY1*C1(K)/(DUMMY**2)
      A11(1,3)=PX1*PZ1*(A1(K)-C1(K))/(DUMMY**2)
      A11(2,1)=A11(1,2)
      A11(2,2)=(PY1/DUMMY)**2*A1(K)+C_DUMMY1**2*B1(K)+(PZ1/DUMMY)**2*S_DUMMY1**2*C1(K)
      A11(2,3)=PY1*PZ1*(A1(K)-C1(K))/(DUMMY**2)
      A11(3,1)=A11(1,3)
      A11(3,2)=A11(2,3)
      A11(3,3)=(PZ1/DUMMY)**2*A1(K)+(PX1**2+PY1**2)*C1(K)/(DUMMY**2)

      DUMMY = SQRT(PX2**2+PY2**2+PZ2**2)

      B11(1,1)=(PX2/DUMMY)**2*A1(K)+S_DUMMY2**2*B1(K)+C_DUMMY2**2*(PZ2/DUMMY)**2*C1(K)
      B11(1,2)=PY2*PX2*A1(K)/(DUMMY**2)-C_DUMMY2*S_DUMMY2*B1(K)+PZ2**2*C_DUMMY2*S_DUMMY2*C1(K)/(DUMMY**2)
      B11(1,3)=PX2*PZ2*(A1(K)-C1(K))/(DUMMY**2)
      B11(2,1)=B11(1,2)
      B11(2,2)=(PY2/DUMMY)**2*A1(K)+C_DUMMY2**2*B1(K)+(PZ2/DUMMY)**2*S_DUMMY2**2*C1(K)
      B11(2,3)=PY2*PZ2*(A1(K)-C1(K))/(DUMMY**2)
      B11(3,1)=B11(1,3)
      B11(3,2)=B11(2,3)
      B11(3,3)=(PZ2/DUMMY)**2*A1(K)+(PX2**2+PY2**2)*C1(K)/(DUMMY**2)


           DO
                  RANDOM=DPRAND()
                   DF_LAMBDA(:)=0.0D0
                   W(:)=0.0D0

                  C_1(1)=(1.0D0-LAMBDA(K))*A11(1,1)+LAMBDA(K)*B11(1,1)
                  C_1(2)=(1.0D0-LAMBDA(K))*A11(1,2)+LAMBDA(K)*B11(1,2)
                  C_1(3)=(1.0D0-LAMBDA(K))*A11(1,3)+LAMBDA(K)*B11(1,3)
                  C_1(4)=(1.0D0-LAMBDA(K))*A11(2,2)+LAMBDA(K)*B11(2,2)
                  C_1(5)=(1.0D0-LAMBDA(K))*A11(2,3)+LAMBDA(K)*B11(2,3)
                  C_1(6)=(1.0D0-LAMBDA(K))*A11(3,3)+LAMBDA(K)*B11(3,3)

 
 
        ! WRITE(*,*) 'CALLING SVERT, MATRIX TO BE INVERTED:'
        ! WRITE(*,*) C_1(:)
        ! WRITE(*,*) ITER
        ! WRITE(*,*) W(:)
          
             CALL SVERT(C_1,3,W)  !INVERTING THE SYMMETRIC MATRIX C1 (OR G, EQ.2.6 IN ECF PAPER)
     

      
             C1INV(1,1)=C_1(1)
             C1INV(1,2)=C_1(2)
             C1INV(1,3)=C_1(3)
             C1INV(2,1)=C1INV(1,2)
             C1INV(2,2)=C_1(4)
             C1INV(2,3)=C_1(5)
             C1INV(3,1)=C1INV(1,3)
             C1INV(3,2)=C1INV(2,3)
             C1INV(3,3)=C_1(6)


        ! WRITE(*,*) 'INVERTED MATRIX:'
        ! WRITE(*,*) C1INV(:,:)
     
            G_INVR=MATMUL(C1INV,R1)
     
             G_INVRT(1,1)=G_INVR(1)
             G_INVRT(1,2)=G_INVR(2)
             G_INVRT(1,3)=G_INVR(3)
     
    
     
             DF_LAMBDA=MATMUL(G_INVRT,MATMUL((1.0D0-LAMBDA(K))**2*A11-LAMBDA(K)**2*B11,G_INVR)) !DERIVATIVE OF THE CONTACT FUNCTION WITH RESPECT TO LAMBDA
             DUMMY=DF_LAMBDA(1)
             DUMMY1=LAMBDA(K)
! WRITE(*,*) "STARTING SUBROUTINE PARAMONOV",ITER, LAMBDA(K), DF_LAMBDA(:)
! WRITE(*,*) G_INVRT(1,:),MATMUL((1.0D0-LAMBDA(K))**2*A11-LAMBDA(K)**2*B11,G_INVR)

             IF(ABS(DF_LAMBDA(1))<CONVG) THEN
                  EXIT
             ELSE
                 ITER=ITER+1

        ! IF THE NEW DERIVATIVE IS LARGER THAN THE OLD ONE, RESET LAMBDA TO ITS PREVIOUS VALUE AND DECREASE STEPSIZE
                IF(ABS(DF_LAMBDA(1))>ABS(DF_LAMBDAOLD)) THEN
                    LAMBDA(K)=LAMBDAOLD
                    DF_LAMBDA(1)=DF_LAMBDAOLD
                    !IF(STEPSIZE>=0.000000001D0)
                    STEPSIZE=STEPSIZE/3.0D0
                END IF

                IF(LAMBDA(K)<0.OR.LAMBDA(K)>1) THEN
!                   WRITE(*,*) "LAMBDA VALUE NOT IN [0,1], RESETTING"
                    STEPSIZE=1.0D0
                    LAMBDA(K)=0.5
                    DF_LAMBDA(1)=0.1D0
                    DUMMY=DF_LAMBDA(1)
                    DUMMY1=LAMBDA(K)
                END IF


                IF(ABS(DF_LAMBDA(1))<=0.1D0) THEN
                    LAMBDA(K)=LAMBDA(K)+DF_LAMBDA(1)*STEPSIZE*RANDOM
                ELSE IF(ABS(DF_LAMBDA(1))>1.0D0) THEN
                    LAMBDA(K)=LAMBDA(K)+(DF_LAMBDA(1)/ABS(DF_LAMBDA(1)))*RANDOM*STEPSIZE*0.2
                ELSE
                    LAMBDA(K)=LAMBDA(K)+DF_LAMBDA(1)*STEPSIZE*0.1*RANDOM
                END IF

 
  LAMBDAOLD=DUMMY1
  DF_LAMBDAOLD=DUMMY

  
             END IF
!           IF(ITER.GT.1000) THEN

!                WRITE(*,*) LAMBDA(K),DF_LAMBDA(1),ITER
!                WRITE(*,*) X1,Y1,Z1
!                WRITE(*,*) X2,Y2,Z2
!                WRITE(*,*) PX1,PY1,PZ1
!                WRITE(*,*) PX2,PY2,PZ2
!                WRITE(*,*) SQRT(PX1**2+PY1**2+PZ1**2)
!                WRITE(*,*) SQRT(PX2**2+PY2**2+PZ2**2)
!           END IF
           END DO
   
    
            S_LAMBDA=LAMBDA(K)*(1.0D0-LAMBDA(K))*MATMUL(R1T,G_INVR)  ! THIS IS THE ELLIPSOID CONTACT FUNCTION
            F1(K)=S_LAMBDA(1)

ECFVALUE=S_LAMBDA(1)

IF(ECFVALUE<1.0D0) THEN
 OVERLAPTEST=.TRUE.
ELSE
 OVERLAPTEST=.FALSE.
END IF

!WRITE(*,*) K
!WRITE(*,*) OVERLAPTEST

END SUBROUTINE OLDECF2



SUBROUTINE GAYBERNE(X,V,EGB,GTEST,STEST)
USE COMMONS
IMPLICIT NONE

! CALCULATION OF THE GAY-BERNE POTENTIAL AND DERIVATIVES
!
! VARIABLES: X1, Y1, Z1, X2, Y2, Z2, PX1, PY1, PZ1, PX2, PY2, PZ2 
! V: 1D ARRAY HOLDING THE DERIVATIVES

LOGICAL   :: GTEST, STEST, OVERLAP, OVERLAPTEST
LOGICAL   :: OVERLAPTESTT
INTEGER   :: I,J1,J2,J3,K1, K2, REALNATOMS
DOUBLE PRECISION :: X(3*NATOMS), V(3*NATOMS), VD(3*NATOMS), X1, Y1, Z1, X2, Y2, Z2, PX1, PY1, PZ1, PX2, PY2, PZ2,&
&    RU1, RU2, R, R2, U1U2, EGB, DPRAND,DERIV, NUMDER(4),DUMMY11, DUMMY22, DUMMY33,&
&    KHI, KHI1, SIGMA, SIGMACUT, FLJ, FLJC, FLJSHIFT, EPSILONNU, EPSILONMU, EPSILON, DUMMY, DUMMY1, &
&    DUMMY2, DUMMY3,DUMMY4, DUMMY5, VDUMMY, RANDOM, &
&    LJ1, LJ2, LJ3, LJ4, LJ6, LJ7, LJ12, LJ13, ECFVALUE, &
&    LJ1C, LJ2C, LJ3C, LJ4C, LJ6C, LJ7C, LJ12C, LJ13C, &
&    LJ1CC, LJ2CC, LJ3CC, LJ4CC, LJ6CC, LJ7CC, LJ12CC, LJ13CC, &
&    RCUT, RCUTU1, RCUTU2, VCUT, &
&    DR_X1, DR_X2, DR_Y1, DR_Y2, DR_Z1, DR_Z2, &
&    DV_PX1, DFLJ_PX1,  DFLJC_PX1, DSIGMA_PX1, DEPSILONNU_PX1, DEPSILONMU_PX1, DEPSILON_PX1, DRU1_PX1, DU1U2_PX1, &
&    DV_PX2, DFLJ_PX2,  DFLJC_PX2, DSIGMA_PX2, DEPSILONNU_PX2, DEPSILONMU_PX2, DEPSILON_PX2, DRU2_PX2, DU1U2_PX2, &
&    DV_PY1, DFLJ_PY1,  DFLJC_PY1, DSIGMA_PY1, DEPSILONNU_PY1, DEPSILONMU_PY1, DEPSILON_PY1, DRU1_PY1, DU1U2_PY1, &
&    DV_PY2, DFLJ_PY2,  DFLJC_PY2, DSIGMA_PY2, DEPSILONNU_PY2, DEPSILONMU_PY2, DEPSILON_PY2, DRU2_PY2, DU1U2_PY2, &
&    DV_PZ1, DFLJ_PZ1,  DFLJC_PZ1, DSIGMA_PZ1, DEPSILONNU_PZ1, DEPSILONMU_PZ1, DEPSILON_PZ1, DRU1_PZ1, DU1U2_PZ1, &
&    DV_PZ2, DFLJ_PZ2,  DFLJC_PZ2, DSIGMA_PZ2, DEPSILONNU_PZ2, DEPSILONMU_PZ2, DEPSILON_PZ2, DRU2_PZ2, DU1U2_PZ2, &
&    DV_X1 , DFLJ_X1 ,  DFLJC_X1 , DSIGMA_X1, DEPSILONNU_X1, DEPSILONMU_X1, DEPSILON_X1, DRU1_X1, DRU2_X1,  &
&    DV_X2 , DFLJ_X2 ,  DFLJC_X2 , DSIGMA_X2, DEPSILONNU_X2, DEPSILONMU_X2, DEPSILON_X2, DRU1_X2, DRU2_X2,  &
&    DV_Y1 , DFLJ_Y1 ,  DFLJC_Y1 , DSIGMA_Y1, DEPSILONNU_Y1, DEPSILONMU_Y1, DEPSILON_Y1, DRU1_Y1, DRU2_Y1,  &
&    DV_Y2 , DFLJ_Y2 ,  DFLJC_Y2 , DSIGMA_Y2, DEPSILONNU_Y2, DEPSILONMU_Y2, DEPSILON_Y2, DRU1_Y2, DRU2_Y2,  &
&    DV_Z1 , DFLJ_Z1 ,  DFLJC_Z1 , DSIGMA_Z1, DEPSILONNU_Z1, DEPSILONMU_Z1, DEPSILON_Z1, DRU1_Z1, DRU2_Z1,  &
&    DV_Z2 , DFLJ_Z2 ,  DFLJC_Z2 , DSIGMA_Z2, DEPSILONNU_Z2, DEPSILONMU_Z2, DEPSILON_Z2, DRU1_Z2, DRU2_Z2,  &
&    KAPPA,KAPPA1     ! SIGMAPAR/SIGMAPERP, ES/EE

DOUBLE PRECISION :: SIGMA0, EPSILON0, MU, NU,  &  ! DEFAULT PARAMETRES OF THE GAY-BERNE POTENTIAL:
&    SIGMAPAR, SIGMAPERP,EE,ES                                ! SIGMAPAR=3.0, MU=2.0, NU=1.0
                                                              ! ES=1.0, EE=0.2: STRENGTH PARAMETERS FOR SIDE-BY-SIDE
                                                              ! AND END-TO-END CONFIGURATIONS

! WRITE(*,*) "GB PARAMETRES:", GBANISOTROPYR, GBWELLDEPTHR, PSIGMA0, PEPSILON0
! WRITE(*,*) 'GAY-BERNE ROUTINE CALLED', OVERLAPTESTT

IF (OVERLAPTESTT) THEN
    OVERLAPTEST = .TRUE.
ELSE
    OVERLAPTEST = .FALSE.
END IF
SIGMA0=PSIGMA0(1)
EPSILON0=PEPSILON0
MU = GBMU
NU = GBNU


SIGMAPAR=GBANISOTROPYR
SIGMAPERP=1.0D0

EE=GBWELLDEPTHR
ES=1.0D0

!A=GBANISOTROPYR/2.0D0
!B=0.5D0

!WRITE(*,*) EE, ES
!WRITE(*,*) SIGMAPAR,SIGMAPERP

! INITIALIZATION OF PARAMETRES

    KHI=(SIGMAPAR**2-SIGMAPERP**2)/(SIGMAPAR**2+SIGMAPERP**2)    ! ANISOTROPY PARAMETER OF ELLIPSOID
    KHI1=(ES**(1.0D0/MU)-EE**(1.0D0/MU))/(ES**(1.0D0/MU)+EE**(1.0D0/MU))   ! ANOTHER ANISOTROPY-LIKE PARAMETER

    KAPPA=SIGMAPAR/SIGMAPERP
    KAPPA1=ES/EE

    RCUT = (GBANISOTROPYR+2.0D0)*SIGMA0!*2.0D0
REALNATOMS=NATOMS/2
I=REALNATOMS

111 CONTINUE

!VT(1:REALNATOMS)=0.0D0
V(1:6*REALNATOMS)=0.0D0
EGB=0.0D0

!WRITE(*,*) 'REALNATOMS', REALNATOMS




DO J1=1,REALNATOMS

    J2=3*J1

    X1 = X(J2-2)
    Y1 = X(J2-1)
    Z1 = X(J2  )

 ! PX, PY, PZ: ANGLE-AXIS COORDINATE SYSTEM, TRYING TO CONVERT THE POTENTIAL TO THESE ONES

    PX1 = X(3*REALNATOMS+J2-2)
    PY1 = X(3*REALNATOMS+J2-1)
    PZ1 = X(3*REALNATOMS+J2  )

!TRY TO NORMALIZE PX, PY, PZ

        DUMMY = SQRT(PX1**2+PY1**2+PZ1**2)
        PX1 = PX1 / DUMMY
        PY1 = PY1 / DUMMY
        PZ1 = PZ1 / DUMMY

! FIRST MAKE PZ1 ALWAYS POSITIVE (DOESN'T MATTER FOR UNIAXIAL ELLIPSOIDS, VEC1 = -VEC1)

!    IF (PZ1.LT.0.0D0) THEN
!        PX1 = -PX1
!        PY1 = -PY1
!        PZ1 = -PZ1
!    END IF


! PUT THE NORMALIZED COORDINATES BACK INTO THE COORDINATE ARRAY
    X(3*REALNATOMS+J2-2) = PX1
    X(3*REALNATOMS+J2-1) = PY1
    X(3*REALNATOMS+J2  ) = PZ1


    DO J3=J1+1,REALNATOMS
  
      VDUMMY=0.0D0
         K2=3*J3
         X2 = X(K2-2)
         Y2 = X(K2-1)
         Z2 = X(K2  )
         PX2 = X(3*REALNATOMS+K2-2)
         PY2 = X(3*REALNATOMS+K2-1)
         PZ2 = X(3*REALNATOMS+K2  )

!! INNER LOOP TO CHECK NUMERICAL DERIVATIVES

!DO I=0,3

!PY2= PY2+0.0001D0
!IF (I==3) PY2 = PY2-0.0004D0


!TRY TO NORMALIZE PX, PY, PZ

!IF(I==0) THEN
        DUMMY = SQRT(PX2**2+PY2**2+PZ2**2)
        PX2 = PX2 / DUMMY
        PY2 = PY2 / DUMMY
        PZ2 = PZ2 / DUMMY
!END IF


! FIRST MAKE PZ2 ALWAYS POSITIVE (DOESN'T MATTER FOR UNIAXIAL ELLIPSOIDS, VEC1 = -VEC1)

!    IF (PZ2.LT.0.0D0) THEN
!        PX2 = -PX2
!        PY2 = -PY2
!        PZ2 = -PZ2
!    END IF

! PUT THE NORMALIZED COORDINATES BACK INTO THE COORDINATE ARRAY
         X(3*REALNATOMS+K2-2) = PX2
         X(3*REALNATOMS+K2-1) = PY2
         X(3*REALNATOMS+K2  ) = PZ2


!         K2=3*K
!         X2=X(K2-2)
!         Y2=X(K2-1)
!         Z2=X(K2  )
!         ALPHA2=X(3*REALNATOMS+K2-2)
!         BETA2=X(3*REALNATOMS+K2-1)
        
! ALPHAR=ATAN((Z2-Z1)/SQRT((X2-X1)**2+(Y2-Y1)**2))          !ALPHAR: ANGLE OF THE R VECTOR AND XY PLANE

        R=SQRT((X2-X1)**2+(Y2-Y1)**2+(Z2-Z1)**2)

! USE CUTOFF
!IF (R>RCUT) THEN
!         VDUMMY = 0.0D0

!         DV_X1=0.0D0
!         DV_X2=0.0D0
!         DV_Y1=0.0D0
!         DV_Y2=0.0D0
!         DV_Z1=0.0D0
!         DV_Z2=0.0D0

!         DV_PX1=0.0D0
!         DV_PY1=0.0D0
!         DV_PZ1=0.0D0

!         DV_PX2=0.0D0
!         DV_PY2=0.0D0
!         DV_PZ2=0.0D0
!WRITE(*,*) R, RCUT
!STOP

!ELSE

! CHECK FOR COLD FUSION
!      IF(ECFVALUE<0.2D0)      THEN   
!         WRITE(*,*) 'SF344> COLD FUSION DETECTED', ECFVALUE 
!         VDUMMY=-1.0D20
!         EGB=-1.0D20
!         V(:)=0.0D0
!         GOTO 999
!      END IF

     
!        IF(OVERLAPTEST) THEN
!                WRITE(*,*) 'GAY-BERNE.F90> OVERLAP DETECTED', ECFVALUE
!        END IF

!      WRITE(*,*) 'VALUE OF THE ECF:', ECFVALUE, OVERLAPTEST

!        WRITE(*,*) 'CHECKING VALUES:', PX2,PY2,PZ2,I


        DR_X1 = (X1-X2)/R
        DR_Y1 = (Y1-Y2)/R
        DR_Z1 = (Z1-Z2)/R
        DR_X2 = -DR_X1
        DR_Y2 = -DR_Y1
        DR_Z2 = -DR_Z1

        DUMMY11 = SQRT(PX1**2+PY1**2+PZ1**2)
        DUMMY22 = SQRT(PX2**2+PY2**2+PZ2**2)
        DUMMY33 = DUMMY11*DUMMY22

        U1U2 = (PX1*PX2 + PY1*PY2 + PZ1*PZ2)/DUMMY33
        RU1 = ((X2-X1)*PX1 + (Y2-Y1)*PY1 + (Z2-Z1)*PZ1)/(R*DUMMY11)
        RU2 = ((X2-X1)*PX2 + (Y2-Y1)*PY2 + (Z2-Z1)*PZ2)/(R*DUMMY22)
!        RCUTU1 =((X2-X1)*PX1 + (Y2-Y1)*PY1 + (Z2-Z1)*PZ1)/(RCUT*DUMMY11)  ! SIGMA DOESN'T DEPEND ON THE ACTUAL DISTANCE
!        RCUTU2 =((X2-X1)*PX2 + (Y2-Y1)*PY2 + (Z2-Z1)*PZ2)/(RCUT*DUMMY22) 

! GET THE VALUE OF SIGMACUT

!        DUMMY1=1.0D0/(1+KHI*U1U2)
!        DUMMY2=1.0D0/(1-KHI*U1U2)
!        DUMMY3=1.0D0-1.0D0/2.0D0*KHI*((RCUTU1+RCUTU2)**2*DUMMY1+(RCUTU1-RCUTU2)**2*DUMMY2)
!        DUMMY4=DUMMY3
!        SIGMACUT= SIGMA0*1.0D0/(SQRT(DUMMY3))
 
!        LJ1CC=1.0D0/(RCUT-SIGMACUT+1.0D0)
!        LJ2CC=LJ1CC**2
!        LJ3CC=LJ2CC*LJ1CC
!        LJ4CC=LJ2CC**2
!        LJ6CC=LJ4CC*LJ2CC
!        LJ7CC=LJ6CC*LJ1CC
!        LJ12CC=LJ6CC**2
!        LJ13CC=LJ12CC*LJ1CC


!((6.0D0*LJ12C - 3.0D0*LJ6C)*((R-SIGMA+1.0D0)/(RCUT-SIGMACUT+1.0D0))**2 - 7.0D0*LJ12C + 4.0D0*LJ6C) 
!               (6.0D0*LJ12C - 3.0D0*LJ6C)*((R-SIGMA+1.0D0)/(RCUT-SIGMACUT+1.0D0))**2 - 7.0D0*LJ12C + 4.0D0*LJ6C

!        DUMMY1=1.0D0/SQRT(1.0D0-KHI**2*(U1U2)**2)
!        EPSILONNU=EPSILON0*DUMMY1
!        DUMMY5=DUMMY1

!        DUMMY1=1.0D0/(1.0D0 + KHI1*U1U2)
!        DUMMY2=1.0D0/(1.0D0 - KHI1*U1U2)
!        EPSILONMU=1.0D0-KHI1/2.0D0*((RCUTU1+RCUTU2)**2*DUMMY1 + (RCUTU1-RCUTU2)**2*DUMMY2)
!        EPSILON=(EPSILONNU**NU)*(EPSILONMU**MU)   

!        LJ1=1.0D0/(RCUT-SIGMA+1.0D0)
!        LJ2=LJ1**2
!        LJ3=LJ2*LJ1
!        LJ4=LJ2**2
!        LJ6=LJ4*LJ2
!        LJ7=LJ6*LJ1
!        LJ12=LJ6**2
!        LJ13=LJ12*LJ1

!        FLJ=LJ12-LJ6

!        VCUT = EPSILON*FLJ

        DUMMY1=1.0D0/(1+KHI*U1U2)
        DUMMY2=1.0D0/(1-KHI*U1U2)
        DUMMY3=1.0D0-1.0D0/2.0D0*KHI*((RU1+RU2)**2*DUMMY1+(RU1-RU2)**2*DUMMY2)
        DUMMY4=DUMMY3
        SIGMA=SIGMA0*1.0D0/(SQRT(DUMMY3))
 
        DUMMY1=1.0D0/SQRT(1.0D0-KHI**2*(U1U2)**2)
        EPSILONNU=EPSILON0*DUMMY1
        DUMMY5=DUMMY1

        DUMMY1=1.0D0/(1.0D0 + KHI1*U1U2)
        DUMMY2=1.0D0/(1.0D0 - KHI1*U1U2)
        EPSILONMU=1.0D0-KHI1/2.0D0*((RU1+RU2)**2*DUMMY1 + (RU1-RU2)**2*DUMMY2)
        EPSILON=(EPSILONNU**NU)*(EPSILONMU**MU)   

! ! LJ-TYPE TERMS FOR R = RCUT
!        LJ1C=1.0D0/(RCUT-SIGMA+1.0D0)
!        LJ2C=LJ1C**2
!        LJ3C=LJ2C*LJ1C
!        LJ4C=LJ2C**2
!        LJ6C=LJ4C*LJ2C
!        LJ7C=LJ6C*LJ1C
!        LJ12C=LJ6C**2
!        LJ13C=LJ12C*LJ1C

! NOW FOR THE ORIGINAL GAY-BERNE POTENTIAL
        LJ1=1.0D0/(R-SIGMA+1.0D0)
        LJ2=LJ1**2
        LJ3=LJ2*LJ1
        LJ4=LJ2**2
        LJ6=LJ4*LJ2
        LJ7=LJ6*LJ1
        LJ12=LJ6**2
        LJ13=LJ12*LJ1

        FLJ=LJ12-LJ6

!   ! STODDARD-FORD TYPE CUTOFF PART OF THE POTENTIAL (SIGMA HAS TO BE CALCULATED FIRST!!! DAMN)
       
!        FLJC =(6.0D0*LJ12C - 3.0D0*LJ6C)*((R-SIGMA+1.0D0)/(RCUT-SIGMA+1.0D0))**2 - 7.0D0*LJ12C + 4.0D0*LJ6C + SIGMA*LJ7CC*(2.0D0*LJ6CC-1.0D0)
!        FLJSHIFT = SIGMA0*LJ7CC*(2.0D0*LJ6CC-1.0D0)

!         WRITE(*,*) 'FLJC:', R, SIGMA, RCUT, SIGMACUT 

        VDUMMY=EPSILON*FLJ      
!        VDUMMY = ((R-SIGMA+1.0D0)/(RCUT-SIGMA+1.0D0))**2 
!        IF(VDUMMY>1.0D5) THEN
!        VDUMMY = -1.0D20
!        V(:) = 1.0D0
!        GOTO 999
!        END IF
 !       WRITE(*,*) 'SF344> ENERGY DEBUG',  DUMMY3

    ! DERIVATIVES FOLLOW 

        DUMMY1 = DUMMY11**3
        DUMMY2 = DUMMY22**3
        DUMMY3 = (PX1*PX2 + PY1*PY2 + PZ1*PZ2)

        DU1U2_PX1 = PX2/DUMMY33 - PX1*DUMMY3/(DUMMY1*DUMMY22)
        DU1U2_PY1 = PY2/DUMMY33 - PY1*DUMMY3/(DUMMY1*DUMMY22)
        DU1U2_PZ1 = PZ2/DUMMY33 - PZ1*DUMMY3/(DUMMY1*DUMMY22)
        DU1U2_PX2 = PX1/DUMMY33 - PX2*DUMMY3/(DUMMY2*DUMMY11)
        DU1U2_PY2 = PY1/DUMMY33 - PY2*DUMMY3/(DUMMY2*DUMMY11)
        DU1U2_PZ2 = PZ1/DUMMY33 - PZ2*DUMMY3/(DUMMY2*DUMMY11)

        DUMMY1 = ((X2-X1)*PX1 + (Y2-Y1)*PY1 + (Z2-Z1)*PZ1)
        DUMMY2 = ((X2-X1)*PX2 + (Y2-Y1)*PY2 + (Z2-Z1)*PZ2)
        R2 = R**2

        DRU1_X1 = (-PX1*R - DR_X1*DUMMY1)/(R2*DUMMY11)
        DRU1_X2 = ( PX1*R - DR_X2*DUMMY1)/(R2*DUMMY11)
        DRU1_Y1 = (-PY1*R - DR_Y1*DUMMY1)/(R2*DUMMY11)
        DRU1_Y2 = ( PY1*R - DR_Y2*DUMMY1)/(R2*DUMMY11)
        DRU1_Z1 = (-PZ1*R - DR_Z1*DUMMY1)/(R2*DUMMY11)
        DRU1_Z2 = ( PZ1*R - DR_Z2*DUMMY1)/(R2*DUMMY11)
        DRU2_X1 = (-PX2*R - DR_X1*DUMMY2)/(R2*DUMMY22)
        DRU2_X2 = ( PX2*R - DR_X2*DUMMY2)/(R2*DUMMY22)
        DRU2_Y1 = (-PY2*R - DR_Y1*DUMMY2)/(R2*DUMMY22)
        DRU2_Y2 = ( PY2*R - DR_Y2*DUMMY2)/(R2*DUMMY22)
        DRU2_Z1 = (-PZ2*R - DR_Z1*DUMMY2)/(R2*DUMMY22)
        DRU2_Z2 = ( PZ2*R - DR_Z2*DUMMY2)/(R2*DUMMY22)


        DRU1_PX1 = (X2-X1)/(R*DUMMY11) - PX1*((X2-X1)*PX1 + (Y2-Y1)*PY1 + (Z2-Z1)*PZ1)/(DUMMY11**3*R)
        DRU1_PY1 = (Y2-Y1)/(R*DUMMY11) - PY1*((X2-X1)*PX1 + (Y2-Y1)*PY1 + (Z2-Z1)*PZ1)/(DUMMY11**3*R)
        DRU1_PZ1 = (Z2-Z1)/(R*DUMMY11) - PZ1*((X2-X1)*PX1 + (Y2-Y1)*PY1 + (Z2-Z1)*PZ1)/(DUMMY11**3*R)
        DRU2_PX2 = (X2-X1)/(R*DUMMY22) - PX2*((X2-X1)*PX2 + (Y2-Y1)*PY2 + (Z2-Z1)*PZ2)/(DUMMY22**3*R)
        DRU2_PY2 = (Y2-Y1)/(R*DUMMY22) - PY2*((X2-X1)*PX2 + (Y2-Y1)*PY2 + (Z2-Z1)*PZ2)/(DUMMY22**3*R)
        DRU2_PZ2 = (Z2-Z1)/(R*DUMMY22) - PZ2*((X2-X1)*PX2 + (Y2-Y1)*PY2 + (Z2-Z1)*PZ2)/(DUMMY22**3*R)


        DUMMY1=1.0D0/(1.0D0 + KHI*U1U2)
        DUMMY2=1.0D0/(1.0D0 - KHI*U1U2)
        DUMMY11 = DUMMY1**2
        DUMMY22 = DUMMY2**2
        
        DSIGMA_PX1=1.0D0/4.0D0*SIGMA0*KHI*1.0D0/SQRT(DUMMY4**3)*((2.0D0*(RU1+RU2)*DRU1_PX1*(1.0D0+KHI*U1U2) &
   &                 -KHI*DU1U2_PX1*(RU1+RU2)**2)*DUMMY11 + &
          & (2.0D0*(RU1-RU2)*DRU1_PX1*(1.0D0-KHI*U1U2)+KHI*DU1U2_PX1*(RU1-RU2)**2)*DUMMY22)
        DSIGMA_PY1=1.0D0/4.0D0*SIGMA0*KHI*1.0D0/SQRT(DUMMY4**3)*((2.0D0*(RU1+RU2)*DRU1_PY1*(1.0D0+KHI*U1U2) &
   &                 -KHI*DU1U2_PY1*(RU1+RU2)**2)*DUMMY11 + &
          & (2.0D0*(RU1-RU2)*DRU1_PY1*(1.0D0-KHI*U1U2)+KHI*DU1U2_PY1*(RU1-RU2)**2)*DUMMY22)
        DSIGMA_PZ1=1.0D0/4.0D0*SIGMA0*KHI*1.0D0/SQRT(DUMMY4**3)*((2.0D0*(RU1+RU2)*DRU1_PZ1*(1.0D0+KHI*U1U2) &
   &                 -KHI*DU1U2_PZ1*(RU1+RU2)**2)*DUMMY11 + &
          & (2.0D0*(RU1-RU2)*DRU1_PZ1*(1.0D0-KHI*U1U2)+KHI*DU1U2_PZ1*(RU1-RU2)**2)*DUMMY22)
 


       DSIGMA_PX2=1.0D0/4.0D0*SIGMA0*KHI*1.0D0/SQRT(DUMMY4**3)*((2.0D0*(RU1+RU2)*DRU2_PX2*(1.0D0+KHI*U1U2) &
   &                               -KHI*DU1U2_PX2*(RU1+RU2)**2)*DUMMY11 + &
      & (2.0D0*(RU1-RU2)*(-DRU2_PX2)*(1.0D0-KHI*U1U2)+KHI*DU1U2_PX2*(RU1-RU2)**2)*DUMMY22)
       DSIGMA_PY2=1.0D0/4.0D0*SIGMA0*KHI*1.0D0/SQRT(DUMMY4**3)*((2.0D0*(RU1+RU2)*DRU2_PY2*(1.0D0+KHI*U1U2) &
   &                               -KHI*DU1U2_PY2*(RU1+RU2)**2)*DUMMY11 + &
      & (2.0D0*(RU1-RU2)*(-DRU2_PY2)*(1.0D0-KHI*U1U2)+KHI*DU1U2_PY2*(RU1-RU2)**2)*DUMMY22)
       DSIGMA_PZ2=1.0D0/4.0D0*SIGMA0*KHI*1.0D0/SQRT(DUMMY4**3)*((2.0D0*(RU1+RU2)*DRU2_PZ2*(1.0D0+KHI*U1U2) &
   &                               -KHI*DU1U2_PZ2*(RU1+RU2)**2)*DUMMY11 + &
      & (2.0D0*(RU1-RU2)*(-DRU2_PZ2)*(1.0D0-KHI*U1U2)+KHI*DU1U2_PZ2*(RU1-RU2)**2)*DUMMY22)


        DSIGMA_X1=1.0D0/4.0D0*SIGMA0*KHI*1.0D0/SQRT(DUMMY4**3)*((2.0D0*(RU1+RU2)*(DRU1_X1+DRU2_X1))*DUMMY1 + &
             & (2.0D0*(RU1-RU2)*(DRU1_X1-DRU2_X1))*DUMMY2)
        DSIGMA_X2=1.0D0/4.0D0*SIGMA0*KHI*1.0D0/SQRT(DUMMY4**3)*((2.0D0*(RU1+RU2)*(DRU1_X2+DRU2_X2))*DUMMY1 + &
             & (2.0D0*(RU1-RU2)*(DRU1_X2-DRU2_X2))*DUMMY2)
        DSIGMA_Y1=1.0D0/4.0D0*SIGMA0*KHI*1.0D0/SQRT(DUMMY4**3)*((2.0D0*(RU1+RU2)*(DRU1_Y1+DRU2_Y1))*DUMMY1 + &
             & (2.0D0*(RU1-RU2)*(DRU1_Y1-DRU2_Y1))*DUMMY2)
        DSIGMA_Y2=1.0D0/4.0D0*SIGMA0*KHI*1.0D0/SQRT(DUMMY4**3)*((2.0D0*(RU1+RU2)*(DRU1_Y2+DRU2_Y2))*DUMMY1 + &
             & (2.0D0*(RU1-RU2)*(DRU1_Y2-DRU2_Y2))*DUMMY2)
        DSIGMA_Z1=1.0D0/4.0D0*SIGMA0*KHI*1.0D0/SQRT(DUMMY4**3)*((2.0D0*(RU1+RU2)*(DRU1_Z1+DRU2_Z1))*DUMMY1 + &
             & (2.0D0*(RU1-RU2)*(DRU1_Z1-DRU2_Z1))*DUMMY2)
        DSIGMA_Z2=1.0D0/4.0D0*SIGMA0*KHI*1.0D0/SQRT(DUMMY4**3)*((2.0D0*(RU1+RU2)*(DRU1_Z2+DRU2_Z2))*DUMMY1 + &
             & (2.0D0*(RU1-RU2)*(DRU1_Z2-DRU2_Z2))*DUMMY2)


         DUMMY1=1.0D0/(1.0D0 + KHI1*U1U2)
         DUMMY2=1.0D0/(1.0D0 - KHI1*U1U2)
         DUMMY11 = DUMMY1**2
         DUMMY22 = DUMMY2**2
    
         DEPSILONMU_PX1=-KHI1/2.0D0*((2.0D0*(RU1+RU2)*DRU1_PX1*(1.0D0+KHI1*U1U2)-KHI1*DU1U2_PX1*(RU1+RU2)**2)*DUMMY11 + &
   &(2.0D0*(RU1-RU2)*DRU1_PX1*(1.0D0-KHI1*U1U2)+KHI1*DU1U2_PX1*(RU1-RU2)**2)*DUMMY22)
         DEPSILONMU_PY1=-KHI1/2.0D0*((2.0D0*(RU1+RU2)*DRU1_PY1*(1.0D0+KHI1*U1U2)-KHI1*DU1U2_PY1*(RU1+RU2)**2)*DUMMY11 + &
   &(2.0D0*(RU1-RU2)*DRU1_PY1*(1.0D0-KHI1*U1U2)+KHI1*DU1U2_PY1*(RU1-RU2)**2)*DUMMY22)
         DEPSILONMU_PZ1=-KHI1/2.0D0*((2.0D0*(RU1+RU2)*DRU1_PZ1*(1.0D0+KHI1*U1U2)-KHI1*DU1U2_PZ1*(RU1+RU2)**2)*DUMMY11 + &
   &(2.0D0*(RU1-RU2)*DRU1_PZ1*(1.0D0-KHI1*U1U2)+KHI1*DU1U2_PZ1*(RU1-RU2)**2)*DUMMY22)

         DEPSILONMU_PX2=-KHI1/2.0D0*((2.0D0*(RU1+RU2)*DRU2_PX2*(1.0D0+KHI1*U1U2)-KHI1*DU1U2_PX2*(RU1+RU2)**2)*DUMMY11 + &
   &(2.0D0*(RU1-RU2)*(-DRU2_PX2)*(1.0D0-KHI1*U1U2)+KHI1*DU1U2_PX2*(RU1-RU2)**2)*DUMMY22)
         DEPSILONMU_PY2=-KHI1/2.0D0*((2.0D0*(RU1+RU2)*DRU2_PY2*(1.0D0+KHI1*U1U2)-KHI1*DU1U2_PY2*(RU1+RU2)**2)*DUMMY11 + &
   &(2.0D0*(RU1-RU2)*(-DRU2_PY2)*(1.0D0-KHI1*U1U2)+KHI1*DU1U2_PY2*(RU1-RU2)**2)*DUMMY22)
         DEPSILONMU_PZ2=-KHI1/2.0D0*((2.0D0*(RU1+RU2)*DRU2_PZ2*(1.0D0+KHI1*U1U2)-KHI1*DU1U2_PZ2*(RU1+RU2)**2)*DUMMY11 + &
   &(2.0D0*(RU1-RU2)*(-DRU2_PZ2)*(1.0D0-KHI1*U1U2)+KHI1*DU1U2_PZ2*(RU1-RU2)**2)*DUMMY22)
 
         DEPSILONMU_X1=-KHI1/2.0D0*(2.0D0*(RU1+RU2)*(DRU1_X1+DRU2_X1)*DUMMY1 + 2.0D0*(RU1-RU2)*(DRU1_X1-DRU2_X1)*DUMMY2)
         DEPSILONMU_X2=-KHI1/2.0D0*(2.0D0*(RU1+RU2)*(DRU1_X2+DRU2_X2)*DUMMY1 + 2.0D0*(RU1-RU2)*(DRU1_X2-DRU2_X2)*DUMMY2)
         DEPSILONMU_Y1=-KHI1/2.0D0*(2.0D0*(RU1+RU2)*(DRU1_Y1+DRU2_Y1)*DUMMY1 + 2.0D0*(RU1-RU2)*(DRU1_Y1-DRU2_Y1)*DUMMY2)
         DEPSILONMU_Y2=-KHI1/2.0D0*(2.0D0*(RU1+RU2)*(DRU1_Y2+DRU2_Y2)*DUMMY1 + 2.0D0*(RU1-RU2)*(DRU1_Y2-DRU2_Y2)*DUMMY2)
         DEPSILONMU_Z1=-KHI1/2.0D0*(2.0D0*(RU1+RU2)*(DRU1_Z1+DRU2_Z1)*DUMMY1 + 2.0D0*(RU1-RU2)*(DRU1_Z1-DRU2_Z1)*DUMMY2)
         DEPSILONMU_Z2=-KHI1/2.0D0*(2.0D0*(RU1+RU2)*(DRU1_Z2+DRU2_Z2)*DUMMY1 + 2.0D0*(RU1-RU2)*(DRU1_Z2-DRU2_Z2)*DUMMY2)
   
         DEPSILONNU_PX1=EPSILON0*KHI**2*DUMMY5**3*U1U2*DU1U2_PX1
         DEPSILONNU_PY1=EPSILON0*KHI**2*DUMMY5**3*U1U2*DU1U2_PY1
         DEPSILONNU_PZ1=EPSILON0*KHI**2*DUMMY5**3*U1U2*DU1U2_PZ1

         DEPSILONNU_PX2=EPSILON0*KHI**2*DUMMY5**3*U1U2*DU1U2_PX2
         DEPSILONNU_PY2=EPSILON0*KHI**2*DUMMY5**3*U1U2*DU1U2_PY2
         DEPSILONNU_PZ2=EPSILON0*KHI**2*DUMMY5**3*U1U2*DU1U2_PZ2

         DEPSILON_PX1=NU*EPSILONNU**(NU-1.0D0)*DEPSILONNU_PX1*EPSILONMU**MU &
   &                       + MU*EPSILONMU**(MU-1.0D0)*DEPSILONMU_PX1*EPSILONNU**NU
         DEPSILON_PY1=NU*EPSILONNU**(NU-1.0D0)*DEPSILONNU_PY1*EPSILONMU**MU & 
   &                       + MU*EPSILONMU**(MU-1.0D0)*DEPSILONMU_PY1*EPSILONNU**NU
         DEPSILON_PZ1=NU*EPSILONNU**(NU-1.0D0)*DEPSILONNU_PZ1*EPSILONMU**MU & 
   &                       + MU*EPSILONMU**(MU-1.0D0)*DEPSILONMU_PZ1*EPSILONNU**NU

         DEPSILON_PX2=NU*EPSILONNU**(NU-1.0D0)*DEPSILONNU_PX2*EPSILONMU**MU &
   &                       + MU*EPSILONMU**(MU-1.0D0)*DEPSILONMU_PX2*EPSILONNU**NU
         DEPSILON_PY2=NU*EPSILONNU**(NU-1.0D0)*DEPSILONNU_PY2*EPSILONMU**MU &
   &                       + MU*EPSILONMU**(MU-1.0D0)*DEPSILONMU_PY2*EPSILONNU**NU
         DEPSILON_PZ2=NU*EPSILONNU**(NU-1.0D0)*DEPSILONNU_PZ2*EPSILONMU**MU &
   &                       + MU*EPSILONMU**(MU-1.0D0)*DEPSILONMU_PZ2*EPSILONNU**NU


         DEPSILON_X1= MU*EPSILONMU**(MU-1.0D0)*DEPSILONMU_X1*EPSILONNU**NU
         DEPSILON_X2= MU*EPSILONMU**(MU-1.0D0)*DEPSILONMU_X2*EPSILONNU**NU
         DEPSILON_Y1= MU*EPSILONMU**(MU-1.0D0)*DEPSILONMU_Y1*EPSILONNU**NU
         DEPSILON_Y2= MU*EPSILONMU**(MU-1.0D0)*DEPSILONMU_Y2*EPSILONNU**NU
         DEPSILON_Z1= MU*EPSILONMU**(MU-1.0D0)*DEPSILONMU_Z1*EPSILONNU**NU
         DEPSILON_Z2= MU*EPSILONMU**(MU-1.0D0)*DEPSILONMU_Z2*EPSILONNU**NU

!        ! FOR THE STODDARD-FORD CUTOFF PART

!        DFLJC_X1 =   2.0D0*(R-SIGMA+1.0D0)*((-1.0D0/R*(X2-X1)-DSIGMA_X1)*(RCUT-SIGMA+1.0D0)+DSIGMA_X1*(R-SIGMA+1.0D0))*&
!                        & (1.0D0/(RCUT-SIGMA+1.0D0))**3*(6.0D0*LJ12C - 3.0D0*LJ6C)-6.0D0*DSIGMA_X1*LJ7C*(2.0D0*LJ6C-1.0D0)
!        DFLJC_X2 =   2.0D0*(R-SIGMA+1.0D0)*(( 1.0D0/R*(X2-X1)-DSIGMA_X2)*(RCUT-SIGMA+1.0D0)+DSIGMA_X2*(R-SIGMA+1.0D0))*&
!                        & (1.0D0/(RCUT-SIGMA+1.0D0))**3*(6.0D0*LJ12C - 3.0D0*LJ6C)-6.0D0*DSIGMA_X2*LJ7C*(2.0D0*LJ6C-1.0D0)
!        DFLJC_Y1 =   2.0D0*(R-SIGMA+1.0D0)*((-1.0D0/R*(Y2-Y1)-DSIGMA_Y1)*(RCUT-SIGMA+1.0D0)+DSIGMA_Y1*(R-SIGMA+1.0D0))*&
!                        & (1.0D0/(RCUT-SIGMA+1.0D0))**3*(6.0D0*LJ12C - 3.0D0*LJ6C)-6.0D0*DSIGMA_Y1*LJ7C*(2.0D0*LJ6C-1.0D0)
!        DFLJC_Y2 =   2.0D0*(R-SIGMA+1.0D0)*(( 1.0D0/R*(Y2-Y1)-DSIGMA_Y2)*(RCUT-SIGMA+1.0D0)+DSIGMA_Y2*(R-SIGMA+1.0D0))*&
!                        & (1.0D0/(RCUT-SIGMA+1.0D0))**3*(6.0D0*LJ12C - 3.0D0*LJ6C)-6.0D0*DSIGMA_Y2*LJ7C*(2.0D0*LJ6C-1.0D0)
!        DFLJC_Z1 =   2.0D0*(R-SIGMA+1.0D0)*((-1.0D0/R*(Z2-Z1)-DSIGMA_Z1)*(RCUT-SIGMA+1.0D0)+DSIGMA_Z1*(R-SIGMA+1.0D0))*&
!                        & (1.0D0/(RCUT-SIGMA+1.0D0))**3*(6.0D0*LJ12C - 3.0D0*LJ6C)-6.0D0*DSIGMA_Z1*LJ7C*(2.0D0*LJ6C-1.0D0)
!        DFLJC_Z2 =   2.0D0*(R-SIGMA+1.0D0)*(( 1.0D0/R*(Z2-Z1)-DSIGMA_Z2)*(RCUT-SIGMA+1.0D0)+DSIGMA_Z2*(R-SIGMA+1.0D0))*&
!                        & (1.0D0/(RCUT-SIGMA+1.0D0))**3*(6.0D0*LJ12C - 3.0D0*LJ6C)-6.0D0*DSIGMA_Z2*LJ7C*(2.0D0*LJ6C-1.0D0)
 
!        DFLJC_PX1 = -2.0D0*(R-SIGMA+1.0D0)*(-DSIGMA_PX1*(RCUT-SIGMA+1.0D0)+DSIGMA_PX1*(R-SIGMA+1))*&
!                        & (1.0D0/(RCUT-SIGMA+1.0D0))**3*(6.0D0*LJ12C-3.0D0*LJ6C)-6.0D0*DSIGMA_PX1*LJ7C*(2.0D0*LJ6C-1.0D0) 
!        DFLJC_PY1 = -2.0D0*(R-SIGMA+1.0D0)*(-DSIGMA_PY1*(RCUT-SIGMA+1.0D0)+DSIGMA_PY1*(R-SIGMA+1))*&
!                        & (1.0D0/(RCUT-SIGMA+1.0D0))**3*(6.0D0*LJ12C-3.0D0*LJ6C)-6.0D0*DSIGMA_PY1*LJ7C*(2.0D0*LJ6C-1.0D0)
!        DFLJC_PZ1 = -2.0D0*(R-SIGMA+1.0D0)*(-DSIGMA_PZ1*(RCUT-SIGMA+1.0D0)+DSIGMA_PZ1*(R-SIGMA+1))*&
!                        & (1.0D0/(RCUT-SIGMA+1.0D0))**3*(6.0D0*LJ12C-3.0D0*LJ6C)-6.0D0*DSIGMA_PZ1*LJ7C*(2.0D0*LJ6C-1.0D0)
!        DFLJC_PX2 = -2.0D0*(R-SIGMA+1.0D0)*(-DSIGMA_PX2*(RCUT-SIGMA+1.0D0)+DSIGMA_PX2*(R-SIGMA+1))*&
!                        & (1.0D0/(RCUT-SIGMA+1.0D0))**3*(6.0D0*LJ12C-3.0D0*LJ6C)-6.0D0*DSIGMA_PX2*LJ7C*(2.0D0*LJ6C-1.0D0)
!        DFLJC_PY2 = -2.0D0*(R-SIGMA+1.0D0)*(-DSIGMA_PY2*(RCUT-SIGMA+1.0D0)+DSIGMA_PY2*(R-SIGMA+1))*&
!                        & (1.0D0/(RCUT-SIGMA+1.0D0))**3*(6.0D0*LJ12C-3.0D0*LJ6C)-6.0D0*DSIGMA_PY2*LJ7C*(2.0D0*LJ6C-1.0D0)
!        DFLJC_PZ2 = -2.0D0*(R-SIGMA+1.0D0)*(-DSIGMA_PZ2*(RCUT-SIGMA+1.0D0)+DSIGMA_PZ2*(R-SIGMA+1))*&
!                        & (1.0D0/(RCUT-SIGMA+1.0D0))**3*(6.0D0*LJ12C-3.0D0*LJ6C)-6.0D0*DSIGMA_PZ2*LJ7C*(2.0D0*LJ6C-1.0D0)

        ! NOW FOR THE ORIGINAL GB POTENTIAL + THE CUTOFF PART

         DFLJ_PX1=6.0D0*DSIGMA_PX1*LJ7*(2.0D0*LJ6-1.0D0)! + DFLJC_PX1
         DFLJ_PY1=6.0D0*DSIGMA_PY1*LJ7*(2.0D0*LJ6-1.0D0)! + DFLJC_PY1
         DFLJ_PZ1=6.0D0*DSIGMA_PZ1*LJ7*(2.0D0*LJ6-1.0D0)! + DFLJC_PZ1
         DFLJ_PX2=6.0D0*DSIGMA_PX2*LJ7*(2.0D0*LJ6-1.0D0)! + DFLJC_PX2
         DFLJ_PY2=6.0D0*DSIGMA_PY2*LJ7*(2.0D0*LJ6-1.0D0)! + DFLJC_PY2
         DFLJ_PZ2=6.0D0*DSIGMA_PZ2*LJ7*(2.0D0*LJ6-1.0D0)! + DFLJC_PZ2

         DFLJ_X1=-6.0D0*LJ7*(2.0D0*LJ6-1.0D0)*(-1.0D0/R*(X2-X1)-DSIGMA_X1)! + DFLJC_X1
         DFLJ_X2=-6.0D0*LJ7*(2.0D0*LJ6-1.0D0)*( 1.0D0/R*(X2-X1)-DSIGMA_X2)! + DFLJC_X2
         DFLJ_Y1=-6.0D0*LJ7*(2.0D0*LJ6-1.0D0)*(-1.0D0/R*(Y2-Y1)-DSIGMA_Y1)! + DFLJC_Y1
         DFLJ_Y2=-6.0D0*LJ7*(2.0D0*LJ6-1.0D0)*( 1.0D0/R*(Y2-Y1)-DSIGMA_Y2)! + DFLJC_Y2
         DFLJ_Z1=-6.0D0*LJ7*(2.0D0*LJ6-1.0D0)*(-1.0D0/R*(Z2-Z1)-DSIGMA_Z1)! + DFLJC_Z1
         DFLJ_Z2=-6.0D0*LJ7*(2.0D0*LJ6-1.0D0)*( 1.0D0/R*(Z2-Z1)-DSIGMA_Z2)! + DFLJC_Z2


         DV_X1=DEPSILON_X1*(FLJ) + EPSILON*DFLJ_X1 
         DV_X2=DEPSILON_X2*(FLJ) + EPSILON*DFLJ_X2 
         DV_Y1=DEPSILON_Y1*(FLJ) + EPSILON*DFLJ_Y1 
         DV_Y2=DEPSILON_Y2*(FLJ) + EPSILON*DFLJ_Y2 
         DV_Z1=DEPSILON_Z1*(FLJ) + EPSILON*DFLJ_Z1 
         DV_Z2=DEPSILON_Z2*(FLJ) + EPSILON*DFLJ_Z2 

         DV_PX1=DEPSILON_PX1*(FLJ) + EPSILON*DFLJ_PX1 
         DV_PY1=DEPSILON_PY1*(FLJ) + EPSILON*DFLJ_PY1 
         DV_PZ1=DEPSILON_PZ1*(FLJ) + EPSILON*DFLJ_PZ1 

         DV_PX2=DEPSILON_PX2*(FLJ) + EPSILON*DFLJ_PX2 
         DV_PY2=DEPSILON_PY2*(FLJ) + EPSILON*DFLJ_PY2 
         DV_PZ2=DEPSILON_PZ2*(FLJ) + EPSILON*DFLJ_PZ2 

!END IF
!        WRITE(*,*) 'SF344: FUNCTION VALUES:'
!        WRITE(*,*) VDUMMY,FLJC,R
!        WRITE(*,*) LJ12-LJ6+((6.0D0*LJ12C - 3.0D0*LJ6C)*((R-SIGMA+1.0D0)/(RCUT-SIGMACUT+1.0D0))**2)- 7.0D0*LJ12C + 4.0D0*LJ6C
!        WRITE(*,*) FLJ +((6.0D0*LJ12C - 3.0D0*LJ6C)*((R-SIGMA+1.0D0)/(RCUT-SIGMACUT+1.0D0))**2- 7.0D0*LJ12C + 4.0D0*LJ6C)
!        WRITE(*,*) SIGMACUT, SIGMA
!        WRITE(*,*) DV_PX1, DV_PX2, DV_PY1, DV_PY2, VDUMMY 
!        STOP

   
         EGB=EGB+VDUMMY
    

!        VT(J) = VT(J) + VDUMMY
!        VT(K) = VT(K) + VDUMMY        ! PAIR POTENTIALS

         V(J2-2)=V(J2-2)+DV_X1
         V(K2-2)=V(K2-2)+DV_X2
         V(J2-1)=V(J2-1)+DV_Y1
         V(K2-1)=V(K2-1)+DV_Y2
         V(J2  )=V(J2  )+DV_Z1
         V(K2  )=V(K2  )+DV_Z2
         V(3*REALNATOMS+J2-2)=V(3*REALNATOMS+J2-2)+DV_PX1
         V(3*REALNATOMS+K2-2)=V(3*REALNATOMS+K2-2)+DV_PX2
         V(3*REALNATOMS+J2-1)=V(3*REALNATOMS+J2-1)+DV_PY1
         V(3*REALNATOMS+K2-1)=V(3*REALNATOMS+K2-1)+DV_PY2
         V(3*REALNATOMS+J2  )=V(3*REALNATOMS+J2  )+DV_PZ1
         V(3*REALNATOMS+K2  )=V(3*REALNATOMS+K2  )+DV_PZ2
!WRITE(*,*) 'SF344> ', V(:) 

!        NUMDER(I+1) = RU2
       
! END DO
!        DERIV = (NUMDER(3)-NUMDER(1))/0.0002
!        WRITE(*,'(A, 3F18.10)') 'SF344> NUMERICAL DERIVATIVES', DERIV, DRU2_PY2, DERIV-DRU2_PY2
!        WRITE(*,*) X(:)
!        STOP
    END DO
END DO
!WRITE(*,*) 'OVERLAPTEST:', OVERLAPTEST
!WRITE(*,*) 'SF344> CHECKING VALUES', R, SIGMA, R-SIGMA
! CHECK FOR OVERLAP
!        WRITE(*,*) X1,Y1,Z1,PX1,PY1,PZ1,X2,Y2,Z2,PX2,PY2,PZ2,GBANISOTROPYR/2.0D0,0.5D0,0.5D0

!        ECFVALUE = 1.0D0

IF(OVERLAPTEST) THEN
!WRITE(*,*) 'IN HERE 1'
DO 
   OVERLAP = .FALSE.

   DO J1 = 1,REALNATOMS
      J2 = 3*J1

      DO K1 = J1+1,REALNATOMS
         K2 = 3*K1
         R = SQRT((X(J2-2)-X(K2-2))**2+(X(J2-1)-X(K2-1))**2+(X(J2)-X(K2))**2)
         IF (R<GBANISOTROPYR) THEN
          DO
           CALL ECF(OVERLAPTEST,ECFVALUE,X(J2-2),X(K2-2),X(J2-1),X(K2-1),X(J2),X(K2),&
                        & X(3*REALNATOMS+J2-2),X(3*REALNATOMS+K2-2),X(3*REALNATOMS+J2-1),X(3*REALNATOMS+K2-1),&
                        & X(3*REALNATOMS+J2),X(3*REALNATOMS+K2),GBANISOTROPYR/2.0D0,0.5D0,0.5D0)
!            WRITE(*,*) 'IN HERE 2', OVERLAPTEST, ECFVALUE
          IF(OVERLAPTEST) THEN
                OVERLAP = .TRUE.
!                RANDOM = DPRAND()
                X(K2-2) = X(K2-2) + 2.0D0*(0.5D0-DPRAND())*0.1
                X(K2-1) = X(K2-1) + 2.0D0*(0.5D0-DPRAND())*0.1
                X(K2  ) = X(K2  ) + 2.0D0*(0.5D0-DPRAND())*0.1
!                X(3*REALNATOMS+K2-2) = X(3*REALNATOMS+K2-2) + 2.0D0*(0.5D0-DPRAND())*0.1
!                X(3*REALNATOMS+K2-2) = X(3*REALNATOMS+K2-2) + 2.0D0*(0.5D0-DPRAND())*0.1
!                X(3*REALNATOMS+K2-2) = X(3*REALNATOMS+K2-2) + 2.0D0*(0.5D0-DPRAND())*0.1
           ELSE
                EXIT
           END IF
          END DO
         END IF
      END DO
   END DO

   IF(.NOT.OVERLAP) THEN
        OVERLAPTEST = .FALSE.
        GOTO 111
   END IF

END DO
END IF
999 CONTINUE











END SUBROUTINE GAYBERNE


SUBROUTINE OLDECF(OVERLAPTEST,ECFVALUE,X1,X2,Y1,Y2,Z1,Z2,ALPHA1,ALPHA2,BETA1,BETA2,A,B)


!CALCULATES THE VALUE OF THE ELLIPSOID CONTACT FUNCTION
USE COMMONS
!USE KEYWORD 
IMPLICIT NONE

DOUBLE PRECISION       ::X1,X2,Y1,Y2,Z1,Z2,ALPHA1,ALPHA2,BETA1,BETA2,A,B
DOUBLE PRECISION  ::LAMBDA,C_ALPHA1,C_BETA1,C_ALPHA2,C_BETA2,S_ALPHA1,S_BETA1,S_ALPHA2,S_BETA2,A1(3,3),B1(3,3),C1(3,3),D1(6),&
&    RVECT(3),RVECTT(1,3),S_LAMBDA(1),ECFVALUE
DOUBLE PRECISION  ::S_LAMBDAMAX(101,1)
INTEGER    ::I,J,K
LOGICAL    ::OVERLAPTEST

S_LAMBDAMAX(:,1)=0.0D0

A=GBANISOTROPYR/2.0D0
B=0.5D0


S_ALPHA1=SIN(ALPHA1)
S_BETA1=SIN(BETA1)
S_ALPHA2=SIN(ALPHA2)
S_BETA2=SIN(BETA2)
C_ALPHA1=COS(ALPHA1)
C_BETA1=COS(BETA1)
C_ALPHA2=COS(ALPHA2)
C_BETA2=COS(BETA2)


A1(1,1) = ((C_ALPHA1*C_BETA1)**2)/A**2 +((S_BETA1)**2)/B**2 +((C_BETA1*S_ALPHA1)**2)/B**2
A1(1,2) = 0
A1(1,3) = 0
A1(2,1) = A1(1,2)
A1(2,2) = ((C_ALPHA1*S_BETA1)**2)/A**2 +((C_BETA1)**2)/B**2 +((S_ALPHA1*S_BETA1)**2)/B**2
A1(2,3) = 0
A1(3,1) = A1(1,3)
A1(3,2) = A1(2,3)
A1(3,3) = (S_ALPHA1)**2/A**2 +(C_ALPHA1)**2/B**2

B1(1,1) = ((C_ALPHA2*C_BETA2)**2)/A**2 +((S_BETA2)**2)/B**2 +((C_BETA2*S_ALPHA2)**2)/B**2
B1(1,2) = 0
B1(1,3) = 0
B1(2,1) = B1(1,2)
B1(2,2) = ((C_ALPHA2*S_BETA2)**2)/A**2 +((C_BETA2)**2)/B**2 +((S_ALPHA2*S_BETA2)**2)/B**2
B1(2,3) = 0
B1(3,1) = B1(1,3)
B1(3,2) = B1(2,3)
B1(3,3) = (S_ALPHA2)**2/A**2 +(C_ALPHA2)**2/B**2

!DIAGONAL MATRIX INVERSION



DO I=1,3
    A1(I,I)=1/A1(I,I)
    B1(I,I)=1/B1(I,I)
END DO



!WRITE(*,*) A1(:,:)
!WRITE(*,*)
!WRITE(*,*) B1(:,:)


RVECT(1)=X2-X1
RVECT(2)=Y2-Y1
RVECT(3)=Z2-Z1
RVECTT(1,1)=RVECT(1)
RVECTT(1,2)=RVECT(2)
RVECTT(1,3)=RVECT(3)

DO K=0,100
    LAMBDA=0.01*K
    
    DO I=1,3
 C1(I,I)=(1-LAMBDA)*A1(I,I)+LAMBDA*B1(I,I)
 C1(I,I)=1/C1(I,I)
    END DO
    
    S_LAMBDA=LAMBDA*(1-LAMBDA)*MATMUL(RVECTT,MATMUL(C1,RVECT))
    S_LAMBDAMAX(K+1,1)=S_LAMBDA(1)
    IF(K>0.AND.S_LAMBDAMAX(K+1,1)<=S_LAMBDAMAX(K,1)) EXIT
END DO



ECFVALUE=S_LAMBDAMAX(K,1)

IF(ECFVALUE<1.1D0) THEN
 OVERLAPTEST=.TRUE.
ELSE
 OVERLAPTEST=.FALSE.
END IF

!WRITE(*,*) K
!WRITE(*,*) OVERLAPTEST

END SUBROUTINE OLDECF


SUBROUTINE OLDECFCHECK(X,GBOVERLAP)
USE COMMONS
!USE KEYWORD
IMPLICIT NONE

DOUBLE PRECISION  :: X(3*NATOMS),X1,X2,Y1,Y2,Z1,Z2,PX1,PX2,PY1,PY2,PZ1,PZ2,ECFVALUE,A,B,C,R
INTEGER    :: J1,J2,J3,J4,REALNATOMS
LOGICAL    :: GBOVERLAP,OVERLAPTEST

IF(GAYBERNET) THEN
        A=(GBANISOTROPYR/2.0D0)**2
        B=(0.5D0)**2
        C=(0.5D0)**2
ELSE IF(PARAMONOVT) THEN
        A=PARAMA1**2
        B=PARAMB1**2
        C=PARAMC1**2
ELSE
 WRITE(*,*) 'ECFCHECK> THIS ROUTINE IS INTENDED TO CHECK OVERLAP BETWEEN ELLIPSOIDS, PLEASE SPECIFY &
&                WHICH SYSTEM YOU WANT TO USE (GAYBERNET OR PARAMONOVT)'
  STOP
            
END IF

REALNATOMS=NATOMS/2

GBOVERLAP=.FALSE.


!WRITE(*,*) A,B

DO J1=1,REALNATOMS
    J2=3*J1
      
          X1=X(J2-2)
          Y1=X(J2-1)
          Z1=X(J2)
          PX1=X(3*REALNATOMS+J2-2)
          PY1=X(3*REALNATOMS+J2-1)
          PZ1=X(3*REALNATOMS+J2  )

  DO J3=J1+1,REALNATOMS
      J4=3*J3
  OVERLAPTEST=.FALSE.
       X2=X(J4-2)
       Y2=X(J4-1)
       Z2=X(J4)
       PX2=X(3*REALNATOMS+J4-2)
       PY2=X(3*REALNATOMS+J4-1)
       PZ2=X(3*REALNATOMS+J4  )

       R=SQRT((X2-X1)**2+(Y2-Y1)**2+(Z2-Z1)**2)
     
      CALL ECF(OVERLAPTEST,ECFVALUE,X1,X2,Y1,Y2,Z1,Z2,PX1,PX2,PY1,PY2,PZ1,PZ2,A,B,C)
!      WRITE(*,*) ECFVALUE, A,B
!      WRITE(*,*) X1,Y1,Z1,X2,Y2,Z2,ALPHA1,BETA1
!      WRITE(*,*) ALPHA2,BETA2
      IF (OVERLAPTEST) THEN
         GBOVERLAP=.TRUE.
!   WRITE(*,*) J1, J3
      END IF
!      WRITE(*,*) OVERLAPTEST
      END DO    
!      WRITE(*,*) GBOVERLAP
!    WRITE(*,*) X(:)
      
END DO


END SUBROUTINE OLDECFCHECK

!
!      ________________________________________________________
!     |                                                        |
!     |       INVERT A SYMMETRIC MATRIX WITHOUT PIVOTING       |
!     |                                                        |
!     |    INPUT:                                              |
!     |                                                        |
!     |        V     --ARRAY PACKED WITH ELEMENTS CONTAINED IN |
!     |                EACH ROW, ON DIAGONAL AND TO RIGHT, OF  |
!     |                COEFFICIENT MATRIX                      |
!     |                (LENGTH AT LEAST N(N+1)/2)              |
!     |                                                        |
!     |        N     --MATRIX DIMENSION                        |
!     |                                                        |
!     |        W     --WORK ARRAY WITH AT LEAST N ELEMENTS     |
!     |                                                        |
!     |    OUTPUT:                                             |
!     |                                                        |
!     |        V     --INVERSE STORED IN THE COMPRESSED MODE   |
!     |                DESCRIBED ABOVE                         |
!     |________________________________________________________|
!
      SUBROUTINE SVERT(V,N,W)
      DOUBLE PRECISION V(*),W(*),S,T
      INTEGER G,H,I,J,K,L,M,N
      H = N
      K = 1
10    IF ( H .EQ. 1 ) GOTO 40
!     --------------------------
!     |*** SAVE PIVOT ENTRY ***|
!     --------------------------
      S = V(K)
      K = K + H
      G = K
      H = H - 1
      M = H
      IF ( S .EQ. 0. ) GOTO 50
      J = 0
20    J = J - M
      M = M - 1
      L = G + M
      T = V(G+J)/S
!     ---------------------------
!     |*** ELIMINATE BY ROWS ***|
!     ---------------------------
      DO 30 I = G,L
30         V(I) = V(I) - T*V(I+J)
      G = L + 1
      IF ( M .GT. 0 ) GOTO 20
      GOTO 10
40    IF ( V(K) .NE. 0. ) GOTO 60
50    WRITE(6,*) 'ERROR: ZERO PIVOT ENCOUNTERED'
      STOP
!     ------------------------------------------
!     |*** SOLVE FOR ROWS OF INVERSE MATRIX ***|
!     ------------------------------------------
60    G = N + N
      DO 150 M = 1,N
           L = ((G-M)*(M-1))/2
           H = L
           K = M
           DO 70 I = M,N
70              W(I) = 0.
           W(M) = 1.
80         IF ( K .EQ. N ) GOTO 100
           T = W(K)/V(K+L)
           J = L
           L = L + N - K
           K = K + 1
           IF ( T .EQ. 0. ) GOTO 80
           DO 90 I = K,N
90              W(I) = W(I) - T*V(I+J)
           GOTO 80
!     -----------------------------------
!     |*** BACK SUBSTITUTION BY ROWS ***|
!     -----------------------------------
100        W(N) = W(N)/V(K+L)
110        IF ( K .EQ. M ) GOTO 130
           J = K
           K = K - 1
           L = L + K - N
           T = W(K)
           DO 120 I = J,N
120             T = T - W(I)*V(I+L)
           W(K) = T/V(K+L)
           GOTO 110
130        DO 140 I = M,N
140             V(I+H) = W(I)
150   CONTINUE
      RETURN
      END SUBROUTINE SVERT

! NUMERIC SECOND DERIVATIVES FOR THE PARAMONOV POTENTIAL


SUBROUTINE PARAMONOVSECDER(OLDX,STEST)
USE COMMONS
!USE KEYWORD
!USE MODHESS
IMPLICIT NONE

DOUBLE PRECISION  :: V(3*NATOMS),EGB,X(3*NATOMS),OLDX(3*NATOMS),VTEMP(2,3*NATOMS),KSI
INTEGER    :: I,J
LOGICAL    :: GTEST,STEST


KSI=0.0001D0
X(:)=OLDX(:)

!WRITE(*,*) 'SF344> IN PARAMONOVSECDER'
DO I=1,3*NATOMS

!       IF ((I.GT.3*NATOMS/2).AND.(MOD(I,3).EQ.0)) THEN
!        VTEMP(:,:) = 0.0D0
!       ELSE
         X(I)=X(I)-KSI
 
         CALL OLDPARAMONOV(X,V,EGB,GTEST,.FALSE.)
 
         VTEMP(1,:)=V(:)
 
         X(I)=X(I)+2.0D0*KSI

         CALL OLDPARAMONOV(X,V,EGB,GTEST,.FALSE.)
 
         VTEMP(2,:)=V(:)

                DO J=1,3*NATOMS
                        HESS(J,I)=(VTEMP(2,J)-VTEMP(1,J))/(2.0D0*KSI)
                END DO
!        END IF
END DO
!WRITE(*,*) 'SF344> EXITING PARAMONOVSECDER'
!WRITE(*,*) 'HESSIAN:'
!WRITE(*,'12F10.3') HESS(:,:)

END SUBROUTINE PARAMONOVSECDER

SUBROUTINE GAYBERNESECDER(OLDX,STEST)
USE COMMONS
!USE KEYWORDS
!USE MODHESS
IMPLICIT NONE

DOUBLE PRECISION  :: V(3*NATOMS),EGB,X(3*NATOMS),OLDX(3*NATOMS),VTEMP(2,3*NATOMS),KSI
INTEGER    :: I,J
LOGICAL    :: GTEST,STEST


KSI=0.000001D0
X(:)=OLDX(:)

!WRITE(*,*) 'SF344> IN GAYBERNESECDER'
DO I=1,3*NATOMS

         X(I)=X(I)-KSI
 
         CALL GAYBERNE(X,V,EGB,GTEST,.FALSE.)
 
         VTEMP(1,:)=V(:)
 
         X(I)=X(I)+2.0D0*KSI

         CALL GAYBERNE(X,V,EGB,GTEST,.FALSE.)
 
         VTEMP(2,:)=V(:)

                DO J=1,3*NATOMS
                        HESS(J,I)=(VTEMP(2,J)-VTEMP(1,J))/(2.0D0*KSI)
                END DO
END DO
!WRITE(*,*) 'SF344> EXITING GAYBERNESECDER'
!WRITE(*,*) 'HESSIAN:'
!WRITE(*,'12F10.3') HESS(:,:)

END SUBROUTINE GAYBERNESECDER


SUBROUTINE TAKESTEPELLIPSOIDS(NP)
! TAKE STEP ROUTINE FOR ANISOTROPIC PARTICLES


USE COMMONS
IMPLICIT NONE

INTEGER                 :: NP, J1, J2, JMAX, J3,K1,K2,K3, JMAX2,J,K,I, REALNATOMS,JP,CLOSESTATOMINDEX
LOGICAL                 :: OVERLAPT, OVERLAPT2, GTEST, DISSOCIATED(NATOMS/2),DISSOC
DOUBLE PRECISION        ::X1,Y1,Z1,X2,Y2,Z2,PX1,PX2,PY1,PY2,PZ1,PZ2,DIST(3*NATOMS/2),R,XMASS,YMASS,ZMASS,DMAX,VMAX,VMAX2
DOUBLE PRECISION        :: VMIN,CMMAX,CMDIST(NATOMS/2),CDIST(NATOMS/2),LOCALSTEP,THETA,PHI,PSI,PI,ECFVALUE,A,B,C,MINDISTANCE

DOUBLE PRECISION        :: EGB, DPRAND,RANDOM,DUMMY, DUMMY1, DUMMY2, DUMMY3,DUMMY4, DUMMY5, RANDOMZ,RVEC(NATOMS/2*(NATOMS/2-1)/2)





!A=PARAMA1**2
!B=PARAMB1**2
!C=PARAMC1**2
IF(GAYBERNET) THEN
        A=(GBANISOTROPYR/2)**2
        B=0.5D0**2
        C=0.5D0**2
ELSE IF(PARAMONOVT) THEN
        A=PARAMA2**2
        B=PARAMB2**2
        C=PARAMC2**2
END IF

PI=ATAN(1.0D0)*4
OVERLAPT=.TRUE.

REALNATOMS=NATOMS/2



    DO J1=1,REALNATOMS
        J2=3*J1
        DUMMY2=COORDS(J2-2,NP)**2+COORDS(J2-1,NP)**2+COORDS(J2,NP)**2
        IF (DUMMY2.GT.RADIUS) THEN
            WRITE(*,'(A,I5,5F20.10)') 'J1,RAD,R**2,X,Y,Z:', J1, RADIUS, DUMMY2, COORDS(J2-2,NP), COORDS(J2-1,NP), COORDS(J2,NP)
            PRINT*, 'INITIAL COORDINATE OUTSIDE CONTAINER -- INCREASE CONTAINER RADIUS'
            STOP
        END IF
    END DO


    DO J1=1,3*NATOMS
        COORDSO(J1,NP)=COORDS(J1,NP)
    END DO

    DO J1=1,NATOMS/2
        VATO(J1,NP)=VAT(J1,NP)
    END DO





!DO J=1,NATOMS
!    WRITE(*,'(A,I2,3F18.10)') "DUMPING COORDINATES", J,COORDS(3*J-2,NP), COORDS(3*J-1,NP), COORDS(3*J,NP)
!END DO

      XMASS=0.0D0
      YMASS=0.0D0
      ZMASS=0.0D0
      DO J1=1,NATOMS/2
         XMASS=XMASS+COORDS(3*(J1-1)+1,NP)
         YMASS=YMASS+COORDS(3*(J1-1)+2,NP)
         ZMASS=ZMASS+COORDS(3*(J1-1)+3,NP)
      ENDDO
      XMASS=XMASS/(NATOMS/2)
      YMASS=YMASS/(NATOMS/2)
      ZMASS=ZMASS/(NATOMS/2)
!
!  FIND THE MOST WEAKLY BOUND ATOM, JMAX, THE SECOND MOST WEAKLY BOUND ATOM, JMAX2,
!  AND THE PAIR ENERGY OF THE MOST TIGHTLY BOUND ATOM, VMIN. AN ANGULAR STEP IS
!  TAKEN FOR JMAX IF ITS PAIR ENERGY IS > ASTEP*VMIN PUTTING THE ATOM AT A RADIUS OF
!  DMAX (OR CMMAX FROM CM OF THE CLUSTER).
!
      DMAX=1.0D0
      VMAX=-1.0D3
      VMAX2=-1.0D3
      VMIN=0.0D0
      CMMAX=1.0D0
      DO J1=1,NATOMS/2
         J2=3*J1
         DIST(J1)= DSQRT( COORDS(J2-2,NP)**2+        COORDS(J2-1,NP)**2+        COORDS(J2,NP)**2)
         CMDIST(J1)=SQRT((COORDS(J2-2,NP)-XMASS)**2+(COORDS(J2-1,NP)-YMASS)**2+(COORDS(J2,NP)-ZMASS)**2)
         IF (CMDIST(J1).GT.CMMAX) CMMAX=CMDIST(J1)
         IF (DIST(J1).GT.DMAX) DMAX=DIST(J1)
         IF (VAT(J1,NP).GT.VMAX) THEN
            VMAX=VAT(J1,NP)
            JMAX=J1
!            WRITE(*,*) "JMAX",JMAX
         ELSE IF ((VAT(J1,NP).LT.VMAX).AND.(VAT(J1,NP).GT.VMAX2)) THEN
            VMAX2=VAT(J1,NP)
            JMAX2=J1
         ENDIF
         IF (VAT(J1,NP).LT.VMIN) VMIN=VAT(J1,NP)
      ENDDO


I=0
OVERLAPT=.FALSE.
DO
I=I+1
OVERLAPT2=.FALSE.
                
            DO J1=1,NATOMS
                J2=3*J1

                LOCALSTEP=STEP(NP)
                IF (J1.GT.NATOMS/2) THEN
                            LOCALSTEP=0.0D0
                        IF (OMOVE(NP)) LOCALSTEP=OSTEP(NP)
                    ELSE IF (J1.LE.NATOMS/2) THEN
                        LOCALSTEP=0.0D0
                        IF (TMOVE(NP)) LOCALSTEP=STEP(NP)
                    END IF



                        IF(VAT(J1,NP).GT.(ASTEP(NP)*VMIN).AND.(J1.EQ.JMAX).AND.(.NOT.NORESET).AND.(J1.LT.NATOMS/2).AND.I==1) THEN
!                                    WRITE(*,*) "ANGULAR MOVE"
                                 THETA=DPRAND()*PI
                                   PHI=DPRAND()*PI*2.0D0
                                   COORDS(J2-2,NP)=XMASS+(CMMAX+1.0D0)*DSIN(THETA)*DCOS(PHI)
                                   COORDS(J2-1,NP)=YMASS+(CMMAX+1.0D0)*DSIN(THETA)*DSIN(PHI)
                                   COORDS(J2,NP)=  ZMASS+(CMMAX+1.0D0)*DCOS(THETA)
                                   DUMMY=COORDS(J2-2,NP)**2+COORDS(J2-1,NP)**2+COORDS(J2,NP)**2
                                   IF (DUMMY.GT.RADIUS) THEN
                                          DUMMY=SQRT(RADIUS*0.99D0/DUMMY)
                                          COORDS(J2-2,NP)=COORDS(J2-2,NP)*DUMMY
                                          COORDS(J2-1,NP)=COORDS(J2-1,NP)*DUMMY
                                          COORDS(J2,NP)=COORDS(J2,NP)*DUMMY
                                END IF
                        ! END OF ANGULAR MOVE BLOCK
                        ELSE IF(I==1) THEN
                            RANDOM=DPRAND()
!                                    PRINT*,'VMIN,VMAX,EFAC=',VMIN,VMAX,EFAC
                            IF ((VMIN-VMAX.EQ.0.0D0).OR.(EFAC.EQ.0.0D0)) THEN
                                    RANDOM=(DPRAND()-0.5D0)*2.0D0
                                    COORDS(J2-2,NP)=COORDS(J2-2,NP)+LOCALSTEP*RANDOM
                                    RANDOM=(DPRAND()-0.5D0)*2.0D0
                                    COORDS(J2-1,NP)=COORDS(J2-1,NP)+LOCALSTEP*RANDOM
                                    RANDOM=(DPRAND()-0.5D0)*2.0D0
                                    COORDS(J2,NP)=COORDS(J2,NP)+LOCALSTEP*RANDOM
                            ENDIF
                        END IF


    
                IF(J1<=REALNATOMS) THEN
                        X1=COORDS(J2-2,NP)
                        Y1=COORDS(J2-1,NP)
                        Z1=COORDS(J2,NP)
                    PX1=COORDS(3*REALNATOMS+J2-2,NP)
                    PY1=COORDS(3*REALNATOMS+J2-1,NP)
                    PZ1=COORDS(3*REALNATOMS+J2  ,NP)

                END IF
        
                DO K=J1+1,NATOMS
                    K2=3*K
                
                IF(K<=REALNATOMS) THEN   
                        X2=COORDS(K2-2,NP)
                        Y2=COORDS(K2-1,NP)
                        Z2=COORDS(K2,NP)
                        PX2=COORDS(3*REALNATOMS+K2-2,NP)
                        PY2=COORDS(3*REALNATOMS+K2-1,NP)
                        PZ2=COORDS(3*REALNATOMS+K2  ,NP)
                       R=SQRT((X2-X1)**2+(Y2-Y1)**2+(Z2-Z1)**2)
                END IF
    
    
                IF(K<=REALNATOMS.AND.J1<=REALNATOMS.AND.R**2<MAX(A,B,C)) THEN
                  CALL ECF(.FALSE.,ECFVALUE,X1,X2,Y1,Y2,Z1,Z2,PX1,PX2,PY1,PY2,PZ1,PZ2,A,B,C)
!                        WRITE(*,*) ECFVALUE                            
                            IF(ECFVALUE<1.0D0)THEN
                              OVERLAPT2=.TRUE.
                    !         WRITE(*,*) "OVERLAPPING ELLIPSOIDS",ECFVALUE,R
                    !         WRITE(*,'(A,3F18.10)') "OLD X,R,SIGMA:",COORDS(K2-2,NP),R,SIGMA
                    !         WRITE(*,'(A,I2,I2)') "J1,K:",J1,K
                              COORDS(K2-2,NP)=COORDS(K2-2,NP)+2*(DPRAND()-0.5D0)*2.0
                              COORDS(K2-1,NP)=COORDS(K2-1,NP)+2*(DPRAND()-0.5D0)*2.0
                              COORDS(K2,NP)=COORDS(K2,NP) +2*(DPRAND()-0.5D0)*2.0
                    !         COORDS(3*REALNATOMS+K2-2,NP)=COORDS(3*REALNATOMS+K2-2,NP)+(DPRAND()-0.5D0)*2.0*ATAN(1.0)*2
                    !         COORDS(3*REALNATOMS+K2-1,NP)=COORDS(3*REALNATOMS+K2-1,NP)+(DPRAND()-0.5D0)*2.0*ATAN(1.0)*2
                            END IF        
                    END IF        
                END DO
            END DO
            
            IF (OVERLAPT2) THEN
                OVERLAPT=.TRUE.
            ELSE
                OVERLAPT=.FALSE.
            END IF


! NOW DETERMINE WHETHER ANY PARTICLE HAS DISSOCIATED
  DO
  DISSOC=.FALSE.
   DO J1=1,NATOMS/2-1
    J2=3*J1
    X1=COORDS(J2-2,NP)
    Y1=COORDS(J2-1,NP)
    Z1=COORDS(J2,NP)
   
    DISSOCIATED(J1)=.TRUE.
    DO K1=J1+1,NATOMS/2
      K2=3*K1
      X2=COORDS(K2-2,NP)
      Y2=COORDS(K2-1,NP)
      Z2=COORDS(K2,NP)
      RVEC(K1)=SQRT((X2-X1)**2+(Y2-Y1)**2+(Z2-Z1)**2)
      IF(RVEC(K1)<PCUTOFF) DISSOCIATED(J1)=.FALSE.
!      WRITE(*,*) 'IN HERE, CUTOFF=',PCUTOFF,J1
    END DO
    IF(DISSOCIATED(J1)) THEN
        ! DETERMINE THE CLOSEST NEIGHBOR
        MINDISTANCE = HUGE(1.0D0)
        DO K1=J1+1,NATOMS/2
         IF(RVEC(K1)<MINDISTANCE) THEN
           MINDISTANCE=RVEC(K1)
           CLOSESTATOMINDEX=K1
         END IF
        END DO
        WRITE(*,*) 'ATOM ', J1, 'IS DISSOCIATED, DISTANCE TO THE CLOSEST ATOM: ', &
        & MINDISTANCE, CLOSESTATOMINDEX
        DISSOC=.TRUE.
        ! MOVE THE DISSOCIATED ATOM CLOSER TO THE CLOSEST ATOM
        COORDS(J2-2,NP)=0.5D0*(COORDS(3*CLOSESTATOMINDEX-2,NP)+COORDS(J2-2,NP))
        COORDS(J2-1,NP)=0.5D0*(COORDS(3*CLOSESTATOMINDEX-1,NP)+COORDS(J2-1,NP))
        COORDS(J2  ,NP)=0.5D0*(COORDS(3*CLOSESTATOMINDEX  ,NP)+COORDS(J2  ,NP))
    END IF
   END DO
  IF(.NOT.DISSOC) EXIT
  END DO
 IF (.NOT.OVERLAPT.AND..NOT.DISSOC) EXIT
END DO                    
!WRITE(*,*) I,"TAKESTEP MOVES"


END SUBROUTINE TAKESTEPELLIPSOIDS

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


SUBROUTINE TAKESTEPPARAM(NP)
! TAKE STEP ROUTINE FOR ANISOTROPIC PARTICLES


USE COMMONS
IMPLICIT NONE

INTEGER                        :: NP, J1, J2, JMAX, J3,K2,K3, JMAX2,J,K,I, REALNATOMS,JP
LOGICAL                        :: OVERLAPT, OVERLAPT2, GTEST

DOUBLE PRECISION        ::X1, Y1, Z1, X2, Y2, Z2, PX1,PY1,PZ1,PX2,PY2,PZ2,DIST(3*NATOMS/2),R,XMASS,YMASS,ZMASS,DMAX,VMAX,VMAX2
DOUBLE PRECISION        :: VMIN,CMMAX,CMDIST(NATOMS/2),CDIST(NATOMS/2),LOCALSTEP,THETA,PHI,PSI,PI,ECFVALUE,A,B,C

DOUBLE PRECISION        :: EGB, DPRAND,RANDOM,DUMMY, DUMMY1, DUMMY2, DUMMY3,DUMMY4, DUMMY5, RANDOMZ




!A=PARAMA1**2
!B=PARAMB1**2
!C=PARAMC1**2

A=PARAMA2**2
B=PARAMB2**2
C=PARAMC2**2


PI=ATAN(1.0D0)*4
OVERLAPT=.TRUE.

REALNATOMS=NATOMS/2



    DO J1=1,REALNATOMS
        J2=3*J1
        DUMMY2=COORDS(J2-2,NP)**2+COORDS(J2-1,NP)**2+COORDS(J2,NP)**2
        IF (DUMMY2.GT.RADIUS) THEN
            WRITE(*,'(A,I5,5F20.10)') 'J1,RAD,R**2,X,Y,Z:', J1, RADIUS, DUMMY2, COORDS(J2-2,NP), COORDS(J2-1,NP), COORDS(J2,NP)
            PRINT*, 'INITIAL COORDINATE OUTSIDE CONTAINER -- INCREASE CONTAINER RADIUS'
            STOP
        END IF
    END DO


    DO J1=1,3*NATOMS
        COORDSO(J1,NP)=COORDS(J1,NP)
    END DO

    DO J1=1,NATOMS/2
        VATO(J1,NP)=VAT(J1,NP)
    END DO





!DO J=1,NATOMS
!    WRITE(*,'(A,I2,3F18.10)') "DUMPING COORDINATES", J,COORDS(3*J-2,NP), COORDS(3*J-1,NP), COORDS(3*J,NP)
!END DO

      XMASS=0.0D0
      YMASS=0.0D0
      ZMASS=0.0D0
      DO J1=1,NATOMS/2
         XMASS=XMASS+COORDS(3*(J1-1)+1,NP)
         YMASS=YMASS+COORDS(3*(J1-1)+2,NP)
         ZMASS=ZMASS+COORDS(3*(J1-1)+3,NP)
      ENDDO
      XMASS=XMASS/(NATOMS/2)
      YMASS=YMASS/(NATOMS/2)
      ZMASS=ZMASS/(NATOMS/2)
!
!  FIND THE MOST WEAKLY BOUND ATOM, JMAX, THE SECOND MOST WEAKLY BOUND ATOM, JMAX2,
!  AND THE PAIR ENERGY OF THE MOST TIGHTLY BOUND ATOM, VMIN. AN ANGULAR STEP IS
!  TAKEN FOR JMAX IF ITS PAIR ENERGY IS > ASTEP*VMIN PUTTING THE ATOM AT A RADIUS OF
!  DMAX (OR CMMAX FROM CM OF THE CLUSTER).
!
      DMAX=1.0D0
      VMAX=-1.0D3
      VMAX2=-1.0D3
      VMIN=0.0D0
      CMMAX=1.0D0
      DO J1=1,NATOMS/2
         J2=3*J1
         DIST(J1)= DSQRT( COORDS(J2-2,NP)**2+        COORDS(J2-1,NP)**2+        COORDS(J2,NP)**2)
         CMDIST(J1)=SQRT((COORDS(J2-2,NP)-XMASS)**2+(COORDS(J2-1,NP)-YMASS)**2+(COORDS(J2,NP)-ZMASS)**2)
         IF (CMDIST(J1).GT.CMMAX) CMMAX=CMDIST(J1)
         IF (DIST(J1).GT.DMAX) DMAX=DIST(J1)
         IF (VAT(J1,NP).GT.VMAX) THEN
            VMAX=VAT(J1,NP)
            JMAX=J1
!            WRITE(*,*) "JMAX",JMAX
         ELSE IF ((VAT(J1,NP).LT.VMAX).AND.(VAT(J1,NP).GT.VMAX2)) THEN
            VMAX2=VAT(J1,NP)
            JMAX2=J1
         ENDIF
         IF (VAT(J1,NP).LT.VMIN) VMIN=VAT(J1,NP)
      ENDDO




I=0
OVERLAPT=.FALSE.
DO
I=I+1
OVERLAPT2=.FALSE.
                
            DO J1=1,NATOMS
                J2=3*J1

                LOCALSTEP=STEP(NP)
                IF (J1.GT.NATOMS/2) THEN
                            LOCALSTEP=0.0D0
                        IF (OMOVE(NP)) LOCALSTEP=OSTEP(NP)
                    ELSE IF (J1.LE.NATOMS/2) THEN
                        LOCALSTEP=0.0D0
                        IF (TMOVE(NP)) LOCALSTEP=STEP(NP)
                    END IF



                        IF(VAT(J1,NP).GT.(ASTEP(NP)*VMIN).AND.(J1.EQ.JMAX).AND.(.NOT.NORESET).AND.(J1.LT.NATOMS/2).AND.I==1) THEN
!                                    WRITE(*,*) "ANGULAR MOVE"
                                 THETA=DPRAND()*PI
                                   PHI=DPRAND()*PI*2.0D0
                                   COORDS(J2-2,NP)=XMASS+(CMMAX+1.0D0)*DSIN(THETA)*DCOS(PHI)
                                   COORDS(J2-1,NP)=YMASS+(CMMAX+1.0D0)*DSIN(THETA)*DSIN(PHI)
                                   COORDS(J2,NP)=  ZMASS+(CMMAX+1.0D0)*DCOS(THETA)
                                   DUMMY=COORDS(J2-2,NP)**2+COORDS(J2-1,NP)**2+COORDS(J2,NP)**2
                                   IF (DUMMY.GT.RADIUS) THEN
                                          DUMMY=SQRT(RADIUS*0.99D0/DUMMY)
                                          COORDS(J2-2,NP)=COORDS(J2-2,NP)*DUMMY
                                          COORDS(J2-1,NP)=COORDS(J2-1,NP)*DUMMY
                                          COORDS(J2,NP)=COORDS(J2,NP)*DUMMY
                                END IF
                        ! END OF ANGULAR MOVE BLOCK
                        ELSE IF(I==1) THEN
                            RANDOM=DPRAND()
!                                    PRINT*,'VMIN,VMAX,EFAC=',VMIN,VMAX,EFAC
                            IF ((VMIN-VMAX.EQ.0.0D0).OR.(EFAC.EQ.0.0D0)) THEN
                                    RANDOM=(DPRAND()-0.5D0)*2.0D0
                                    COORDS(J2-2,NP)=COORDS(J2-2,NP)+LOCALSTEP*RANDOM
                                    RANDOM=(DPRAND()-0.5D0)*2.0D0
                                    COORDS(J2-1,NP)=COORDS(J2-1,NP)+LOCALSTEP*RANDOM
                                    RANDOM=(DPRAND()-0.5D0)*2.0D0
                                    COORDS(J2,NP)=COORDS(J2,NP)+LOCALSTEP*RANDOM
                            ENDIF
                        END IF


    
                IF(J1<=REALNATOMS) THEN
                        X1=COORDS(J2-2,NP)
                        Y1=COORDS(J2-1,NP)
                        Z1=COORDS(J2,NP)
                    PX1=COORDS(3*REALNATOMS+J2-2,NP)
                    PY1=COORDS(3*REALNATOMS+J2-1,NP)
                    PZ1=COORDS(3*REALNATOMS+J2  ,NP)
                END IF
        
                DO K=J1+1,NATOMS
                    K2=3*K
                
                    IF(K<=REALNATOMS) THEN    
                            X2=COORDS(K2-2,NP)
                        Y2=COORDS(K2-1,NP)
                        Z2=COORDS(K2,NP)
                        PX2=COORDS(3*REALNATOMS+K2-2,NP)
                        PY2=COORDS(3*REALNATOMS+K2-1,NP)
                        PZ2=COORDS(3*REALNATOMS+K2  ,NP)

                            R=DSQRT((X2-X1)**2+(Y2-Y1)**2+(Z2-Z1)**2)
                    END IF
                    
                    
                    IF(K<=REALNATOMS.AND.J1<=REALNATOMS.AND.R**2<MAX(A,B,C)) THEN
                                            
                            CALL ECF(.FALSE.,ECFVALUE,X1,X2,Y1,Y2,Z1,Z2,PX1,PX2,PY1,PY2,PZ1,PZ2,A,B,C)
!                        WRITE(*,*) ECFVALUE                            
                            IF(ECFVALUE<1.0D0)THEN
                                            OVERLAPT2=.TRUE.
                    !                            WRITE(*,*) "OVERLAPPING ELLIPSOIDS",ECFVALUE,R
                    !                                WRITE(*,'(A,3F18.10)') "OLD X,R,SIGMA:",COORDS(K2-2,NP),R,SIGMA
                    !                                WRITE(*,'(A,I2,I2)') "J1,K:",J1,K
                            
                                                COORDS(K2-2,NP)=COORDS(K2-2,NP)+2*(DPRAND()-0.5D0)*2.0
                                                   COORDS(K2-1,NP)=COORDS(K2-1,NP)+2*(DPRAND()-0.5D0)*2.0
                                                COORDS(K2,NP)=COORDS(K2,NP) +2*(DPRAND()-0.5D0)*2.0
                    !                            COORDS(3*REALNATOMS+K2-2,NP)=COORDS(3*REALNATOMS+K2-2,NP)+(DPRAND()-0.5D0)*2.0*ATAN(1.0)*2
                        !                            COORDS(3*REALNATOMS+K2-1,NP)=COORDS(3*REALNATOMS+K2-1,NP)+(DPRAND()-0.5D0)*2.0*ATAN(1.0)*2
                                            
                                            

                    !                            WRITE(*,'(A)') "MOVING ORIENTATION"
                                                                            
                                END IF        
                    END IF        
                END DO
            END DO
            
            IF (OVERLAPT2) THEN
                OVERLAPT=.TRUE.
            ELSE
                OVERLAPT=.FALSE.
            END IF
            IF (.NOT.OVERLAPT) EXIT

END DO                    
                    
!WRITE(*,*) I,"TAKESTEP MOVES"


END SUBROUTINE TAKESTEPPARAM


SUBROUTINE TAKESTEPGB(NP)
! TAKE STEP ROUTINE FOR GAY-BERNE ELLIPSOIDS


USE COMMONS
IMPLICIT NONE

INTEGER                        :: NP, J1, J2, JMAX, J3,K2,K3, JMAX2,J,K,I, REALNATOMS,JP
LOGICAL                        :: OVERLAPT, OVERLAPT2, GTEST

DOUBLE PRECISION        ::X1, Y1, Z1, X2, Y2, Z2, PX1,PX2,PY1,PY2,PZ1,PZ2, ALPHAR, &
&                       BETAR,DIST(3*NATOMS/2),XMASS,YMASS,ZMASS,DMAX,VMAX,VMAX2,&
&                        S_ALPHA1,S_ALPHA2,S_ALPHAR,S_BETA12,S_BETA1R,S_BETA2R,C_ALPHA1,C_ALPHA2,&
&                       C_ALPHAR,C_BETA12,C_BETA1R,C_BETA2R,&
&                        VMIN,CMMAX,CMDIST(NATOMS/2),CDIST(NATOMS/2),LOCALSTEP,THETA,PHI,PSI,PI,ECFVALUE,A,B,C

DOUBLE PRECISION        :: R, DPRAND,RANDOM,&
&                         KHI, KHI1, DUMMY, DUMMY1, DUMMY2, DUMMY3,DUMMY4, DUMMY5, RANDOMZ


A=(GBANISOTROPYR/2.0D0)**2
B=(0.5D0)**2
C=(0.5D0)**2

!WRITE(*,*) "GBANISOTROPYR",A

PI=ATAN(1.0D0)*4
OVERLAPT=.TRUE.

REALNATOMS=NATOMS/2



    DO J1=1,REALNATOMS
        J2=3*J1
        DUMMY2=COORDS(J2-2,NP)**2+COORDS(J2-1,NP)**2+COORDS(J2,NP)**2
        IF (DUMMY2.GT.RADIUS) THEN
            WRITE(*,'(A,I5,5F20.10)') 'J1,RAD,R**2,X,Y,Z:', J1, RADIUS, DUMMY2, COORDS(J2-2,NP), COORDS(J2-1,NP), COORDS(J2,NP)
            PRINT*, 'INITIAL COORDINATE OUTSIDE CONTAINER -- INCREASE CONTAINER RADIUS'
            STOP
        END IF
    END DO


    DO J1=1,3*NATOMS
        COORDSO(J1,NP)=COORDS(J1,NP)
    END DO

    DO J1=1,NATOMS/2
        VATO(J1,NP)=VAT(J1,NP)
    END DO






      XMASS=0.0D0
      YMASS=0.0D0
      ZMASS=0.0D0
      DO J1=1,NATOMS/2
         XMASS=XMASS+COORDS(3*(J1-1)+1,NP)
         YMASS=YMASS+COORDS(3*(J1-1)+2,NP)
         ZMASS=ZMASS+COORDS(3*(J1-1)+3,NP)
      ENDDO
      XMASS=XMASS/(NATOMS/2)
      YMASS=YMASS/(NATOMS/2)
      ZMASS=ZMASS/(NATOMS/2)
!
!  FIND THE MOST WEAKLY BOUND ATOM, JMAX, THE SECOND MOST WEAKLY BOUND ATOM, JMAX2,
!  AND THE PAIR ENERGY OF THE MOST TIGHTLY BOUND ATOM, VMIN. AN ANGULAR STEP IS
!  TAKEN FOR JMAX IF ITS PAIR ENERGY IS > ASTEP*VMIN PUTTING THE ATOM AT A RADIUS OF
!  DMAX (OR CMMAX FROM CM OF THE CLUSTER).
!
      DMAX=1.0D0
      VMAX=-1.0D3
      VMAX2=-1.0D3
      VMIN=0.0D0
      CMMAX=1.0D0
      DO J1=1,NATOMS/2
         J2=3*J1
         DIST(J1)= DSQRT( COORDS(J2-2,NP)**2+        COORDS(J2-1,NP)**2+        COORDS(J2,NP)**2)
         CMDIST(J1)=SQRT((COORDS(J2-2,NP)-XMASS)**2+(COORDS(J2-1,NP)-YMASS)**2+(COORDS(J2,NP)-ZMASS)**2)
         IF (CMDIST(J1).GT.CMMAX) CMMAX=CMDIST(J1)
         IF (DIST(J1).GT.DMAX) DMAX=DIST(J1)
         IF (VAT(J1,NP).GT.VMAX) THEN
            VMAX=VAT(J1,NP)
            JMAX=J1
!            WRITE(*,*) "JMAX",JMAX
         ELSE IF ((VAT(J1,NP).LT.VMAX).AND.(VAT(J1,NP).GT.VMAX2)) THEN
            VMAX2=VAT(J1,NP)
            JMAX2=J1
         ENDIF
         IF (VAT(J1,NP).LT.VMIN) VMIN=VAT(J1,NP)
      ENDDO







I=0
OVERLAPT=.FALSE.
DO
I=I+1
OVERLAPT2=.FALSE.
                
            DO J1=1,NATOMS
                J2=3*J1

                LOCALSTEP=STEP(NP)
                IF (J1.GT.NATOMS/2) THEN
                            LOCALSTEP=0.0D0
                        IF (OMOVE(NP)) LOCALSTEP=OSTEP(NP)
                    ELSE IF (J1.LE.NATOMS/2) THEN
                        LOCALSTEP=0.0D0
                        IF (TMOVE(NP)) LOCALSTEP=STEP(NP)
                    END IF



                        IF(VAT(J1,NP).GT.(ASTEP(NP)*VMIN).AND.(J1.EQ.JMAX).AND.(.NOT.NORESET).AND.(J1.LT.NATOMS/2).AND.I==1) THEN
!                                    WRITE(*,*) "ANGULAR MOVE"
                                 THETA=DPRAND()*PI
                                   PHI=DPRAND()*PI*2.0D0
                                   COORDS(J2-2,NP)=XMASS+(CMMAX+1.0D0)*DSIN(THETA)*DCOS(PHI)
                                   COORDS(J2-1,NP)=YMASS+(CMMAX+1.0D0)*DSIN(THETA)*DSIN(PHI)
                                   COORDS(J2,NP)=  ZMASS+(CMMAX+1.0D0)*DCOS(THETA)
                                   DUMMY=COORDS(J2-2,NP)**2+COORDS(J2-1,NP)**2+COORDS(J2,NP)**2
                                   IF (DUMMY.GT.RADIUS) THEN
                                          DUMMY=SQRT(RADIUS*0.99D0/DUMMY)
                                          COORDS(J2-2,NP)=COORDS(J2-2,NP)*DUMMY
                                          COORDS(J2-1,NP)=COORDS(J2-1,NP)*DUMMY
                                          COORDS(J2,NP)=COORDS(J2,NP)*DUMMY
                                END IF
                        ! END OF ANGULAR MOVE BLOCK
                        ELSE IF(I==1) THEN
                            RANDOM=DPRAND()
!                                    PRINT*,'VMIN,VMAX,EFAC=',VMIN,VMAX,EFAC
                            IF ((VMIN-VMAX.EQ.0.0D0).OR.(EFAC.EQ.0.0D0)) THEN
                                    RANDOM=(DPRAND()-0.5D0)*2.0D0
                                    COORDS(J2-2,NP)=COORDS(J2-2,NP)+LOCALSTEP*RANDOM
                                    RANDOM=(DPRAND()-0.5D0)*2.0D0
                                    COORDS(J2-1,NP)=COORDS(J2-1,NP)+LOCALSTEP*RANDOM
                                    RANDOM=(DPRAND()-0.5D0)*2.0D0
                                    COORDS(J2,NP)=COORDS(J2,NP)+LOCALSTEP*RANDOM
                            ENDIF
                        END IF


    
                IF(J1<=REALNATOMS) THEN
                        X1=COORDS(J2-2,NP)
                        Y1=COORDS(J2-1,NP)
                        Z1=COORDS(J2,NP)
                        PX1=COORDS(3*REALNATOMS+J2-2,NP)
                        PY1=COORDS(3*REALNATOMS+J2-1,NP)
                    PZ1=COORDS(3*REALNATOMS+J2  ,NP)
                END IF
        
                DO K=J1+1,NATOMS
                    K2=3*K
                
                    IF(K<=REALNATOMS) THEN    
                            X2=COORDS(K2-2,NP)
                        Y2=COORDS(K2-1,NP)
                        Z2=COORDS(K2,NP)
                        PX2=COORDS(3*REALNATOMS+K2-2,NP)
                        PY2=COORDS(3*REALNATOMS+K2-1,NP)
                        PZ2=COORDS(3*REALNATOMS+K2  ,NP)
                            R=DSQRT((X2-X1)**2+(Y2-Y1)**2+(Z2-Z1)**2)
                    END IF
                    
                    
                                    
                    IF(K<=REALNATOMS.AND.J1<=REALNATOMS.AND.R<=MAX(GBANISOTROPYR,1.0D0)) THEN
                        
                            
                            CALL ECF(.FALSE.,ECFVALUE,X1,X2,Y1,Y2,Z1,Z2,PX1,PX2,PY1,PY2,PZ1,PZ2,A,B,C)
                                        
                            IF(ECFVALUE<1.0D0)THEN
                                            OVERLAPT2=.TRUE.
                    !                            WRITE(*,*) "OVERLAPPING ELLIPSOIDS",ECFVALUE,R
                    !                                WRITE(*,'(A,3F18.10)') "OLD X,R,SIGMA:",COORDS(K2-2,NP),R,SIGMA
                    !                                WRITE(*,'(A,I2,I2)') "J1,K:",J1,K
                            
                                                COORDS(K2-2,NP)=COORDS(K2-2,NP)+2*(DPRAND()-0.5D0)*2.0
                                                   COORDS(K2-1,NP)=COORDS(K2-1,NP)+2*(DPRAND()-0.5D0)*2.0
                                                COORDS(K2,NP)=COORDS(K2,NP) +2*(DPRAND()-0.5D0)*2.0
                    !                            COORDS(3*REALNATOMS+K2-2,NP)=COORDS(3*REALNATOMS+K2-2,NP)+(DPRAND()-0.5D0)*2.0*ATAN(1.0)*2
                        !                            COORDS(3*REALNATOMS+K2-1,NP)=COORDS(3*REALNATOMS+K2-1,NP)+(DPRAND()-0.5D0)*2.0*ATAN(1.0)*2
                                            
                                            

                    !                            WRITE(*,'(A)') "MOVING ORIENTATION"
                                                                            
                                END IF        
                    END IF        
                END DO
            END DO
            
            IF (OVERLAPT2) THEN
                OVERLAPT=.TRUE.
            ELSE
                OVERLAPT=.FALSE.
            END IF
            IF (.NOT.OVERLAPT) EXIT




END DO                    
                    
!WRITE(*,*) I



END SUBROUTINE TAKESTEPGB

!
!	DWAIPAYAN'S STUFF IN HERE
!
!     _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _

      SUBROUTINE GBINIT(X)

      USE COMMONS, ONLY : NATOMS, GBKAPPA, GBKAPPRM, GBMU, GBNU, GBCHI, GBCHIPRM, SIGNOT, EPSNOT
      IMPLICIT NONE

      INTEGER          :: REALNATOMS, OFFSET, I, J
      DOUBLE PRECISION :: X(3*NATOMS), V(3*NATOMS), EGB
      DOUBLE PRECISION :: THETA, PHI, CTI, STI, CPI, SPI, PI
      LOGICAL          :: GTEST
        
!     CALCULATE FROM INPUT PARAMETERS

      REALNATOMS = NATOMS / 2
      OFFSET     = 3 * REALNATOMS

      GBCHI    = (GBKAPPA ** 2 - 1.D0) / (GBKAPPA ** 2 + 1.D0)
      GBCHIPRM = (GBKAPPRM**(1.D0/GBMU)-1.D0) / (GBKAPPRM**(1.D0/GBMU)+1.D0)
      PI       = 4.D0*DATAN(1.D0)

      DO I = 1, REALNATOMS

         J      = OFFSET + 3*I 
         THETA  = DACOS(X(J))

         IF(THETA == 0.0D0 .OR. THETA == PI) THEN 
            PHI = 0.D0
         ELSE
            CPI    = X(J-2) / DSIN(THETA) 
            SPI    = X(J-1) / DSIN(THETA)
            PHI    = DATAN2(SPI,CPI)
            IF (PHI .LT. 0.D0) PHI = PHI + 2.D0*PI
         ENDIF
         X(J-2) = THETA 
         X(J-1) = PHI
         X(J)   = 0.D0

      ENDDO

      END SUBROUTINE GBINIT

!     _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
       
      SUBROUTINE GBGTP(X, V, EGB, GTEST) 

      USE COMMONS, ONLY : NATOMS, GBKAPPA, GBKAPPRM, GBMU, GBNU, GBCHI, GBCHIPRM, SIGNOT, EPSNOT

      IMPLICIT NONE

      INTEGER          :: REALNATOMS, OFFSET
      INTEGER          I, J, J1, J2, J3, J4, J5, J6, K
      DOUBLE PRECISION X(3*NATOMS), V(3*NATOMS)
      DOUBLE PRECISION GX(NATOMS/2), GY(NATOMS/2), GZ(NATOMS/2), EGB
      DOUBLE PRECISION THETA(NATOMS/2), PHI(NATOMS/2)
      DOUBLE PRECISION SCSIG, SCSIG3 
      DOUBLE PRECISION EPS1, EPS2, EPS
      DOUBLE PRECISION FCT1, FCT2, FCT3, FCT4, FCT5, FCT6, FCT7, FCT8
      DOUBLE PRECISION FCT9, FCT10, FCT11, FCT12, FCT13, FCT14, FCT15
      DOUBLE PRECISION FCT16, FCT17, FCT18, FCT19, FCT20
      DOUBLE PRECISION FCT3P4, FCT3M4, FCT7P8, FCT7M8
      DOUBLE PRECISION ALP, BET, GAM, APB, AMB
      DOUBLE PRECISION R, INVR, NR(3), SCR 
      DOUBLE PRECISION SRM1, SRM2, SRM6, SRM7, SRM12, SRM13, SR12M6   
      DOUBLE PRECISION DX, DY, DZ, R2, VR, VA, VB, VG
      DOUBLE PRECISION FXIJ, FYIJ, FZIJ, FIJN, FIJEI, FIJEJ
      DOUBLE PRECISION CTI, STI, CPI, SPI, GI(3), GJ(3), T1, P1
      LOGICAL          GTEST

      REALNATOMS = NATOMS/2
      OFFSET     = 3 * REALNATOMS

!     INITIALIZE 

      EGB               = 0.D0
      V(1:6*REALNATOMS) = 0.D0
      GX(1:REALNATOMS)  = 0.D0
      GY(1:REALNATOMS)  = 0.D0
      GZ(1:REALNATOMS)  = 0.D0

      IF (GTEST) THEN

         DO I = 1, REALNATOMS

            J             = 3*I
            THETA(I)      = X(OFFSET+J-2)
            PHI(I)        = X(OFFSET+J-1)
            T1            = THETA(I)
            P1            = PHI(I)
            X(OFFSET+J-2) = DCOS(P1)*DSIN(T1)
            X(OFFSET+J-1) = DSIN(P1)*DSIN(T1)
            X(OFFSET+J)   = DCOS(T1)

         ENDDO         

         DO J1 = 1, REALNATOMS-1

            J3 = 3*J1 
            J5 = OFFSET + J3
 
            DO J2 = J1 + 1, REALNATOMS 

               J4      = 3*J2
               J6      = OFFSET + J4 

               DX      = X(J3-2) - X(J4-2)
               DY      = X(J3-1) - X(J4-1)
               DZ      = X(J3) - X(J4)
               R2      = DX*DX + DY*DY + DZ*DZ
               R       = DSQRT(R2)
               SCR     = R/SIGNOT
               INVR    = 1.D0/R
               R2      = 1.D0/R2
               
!     NORMILIZE THE SEPRATION VECTOR 

               NR(1)   = DX*INVR
               NR(2)   = DY*INVR
               NR(3)   = DZ*INVR

!     CALCULATE $\ALPHA$, $\BETA$ AND $\GAMMA$

               ALP = NR(1)*X(J5-2)+NR(2)*X(J5-1)+NR(3)*X(J5)
               BET = NR(1)*X(J6-2)+NR(2)*X(J6-1)+NR(3)*X(J6)
               GAM = X(J5-2)*X(J6-2)+X(J5-1)*X(J6-1)+X(J5)*X(J6)

!     CALCULATE USEFUL QUANTITIES

               APB    = ALP+BET
               AMB    = ALP-BET

               FCT1   = 1.D0/(1.D0+GBCHI*GAM)
               FCT2   = 1.D0/(1.D0-GBCHI*GAM)
               FCT3   = APB*FCT1
               FCT4   = AMB*FCT2
               FCT3P4 = FCT3+FCT4
               FCT3M4 = FCT3-FCT4

               FCT5   = 1.D0/(1.D0+GBCHIPRM*GAM)
               FCT6   = 1.D0/(1.D0-GBCHIPRM*GAM)
               FCT7   = (ALP+BET)*FCT5
               FCT8   = (ALP-BET)*FCT6
               FCT7P8 = FCT7+FCT8
               FCT7M8 = FCT7-FCT8

!     CALCULATE $\EPSILON$

               EPS1   = DSQRT(FCT1*FCT2)
               EPS2   = 1.D0-0.5D0*GBCHIPRM*(APB*FCT7+AMB*FCT8)
               EPS    = EPSNOT*EPS1**GBNU*EPS2**GBMU

!     CALCULATE $(\SIGMA/\SIGMA_{0})^3$

               SCSIG  = 1.D0/DSQRT(1.D0-0.5D0*GBCHI*(APB*FCT3+AMB*FCT4))
               SCSIG3 = SCSIG*SCSIG*SCSIG 

!     CALCULATE DEL(V)/DEL(R)

               SRM1   = 1.D0/(SCR-SCSIG+1.D0)  
               SRM2   = SRM1*SRM1              
               SRM6   = SRM2*SRM2*SRM2
               SRM7   = SRM6*SRM1 
               SRM12  = SRM6*SRM6
               SRM13  = SRM12*SRM1
               SR12M6 = SRM12-SRM6
               FCT9   = 2.D0*SRM13-SRM7
               VR     = -(24.D0/SIGNOT)*EPS*FCT9 

!     CALCULATE ENERGY

               EGB = EGB + EPS*SR12M6
               
!     CALCULATE DEL(V)/DEL(\ALPHA) AND DEL(V)/DEL(\BETA) 

               FCT10  = 24.D0*GBMU*GBCHIPRM*EPS/(SIGNOT*EPS2)
               FCT11  = 4.D0*EPS*GBMU*GBCHIPRM/EPS2
               FCT12  = 12.D0*EPS*GBCHI*SCSIG3

               VA     = -FCT11*SR12M6*FCT7P8+FCT12*FCT9*FCT3P4
               VB     = -FCT11*SR12M6*FCT7M8+FCT12*FCT9*FCT3M4

!     CALCULATE DEL(V)/DEL(\GAMMA)  

               FCT13  = EPS1*EPS1*GBCHI*GBCHI*GAM*GBNU
               FCT14  = 0.5D0*GBMU*GBCHIPRM*GBCHIPRM*(FCT7*FCT7-FCT8*FCT8)/EPS2
!     %                  /EPS2
               FCT15  = FCT13 + FCT14
               FCT16  = FCT3*FCT3 - FCT4*FCT4
               FCT17  = 4.D0*EPS*FCT15 
               FCT18  = 6.D0*EPS*GBCHI*GBCHI*SCSIG3*FCT16
               FCT19  = FCT17*SR12M6
               FCT20  = FCT18*FCT9

               VG     = FCT19 - FCT20

               DO K = 1, 3

                  GI(K) = VA*NR(K) + VG*X(J6-3+K)
                  GJ(K) = VB*NR(K) + VG*X(J5-3+K)

               ENDDO           

!     CALCULATE CONTRIBUTION TO FORCES 

               FIJN   = VR - (VA*ALP+VB*BET)*INVR  
               FIJEI  = VA*INVR
               FIJEJ  = VB*INVR

               FXIJ = FIJN*NR(1)+FIJEI*X(J5-2)+FIJEJ*X(J6-2)
               FYIJ = FIJN*NR(2)+FIJEI*X(J5-1)+FIJEJ*X(J6-1)
               FZIJ = FIJN*NR(3)+FIJEI*X(J5)+FIJEJ*X(J6) 

               V(J3-2) = V(J3-2) + FXIJ
               V(J3-1) = V(J3-1) + FYIJ
               V(J3)   = V(J3) + FZIJ

               V(J4-2) = V(J4-2) - FXIJ
               V(J4-1) = V(J4-1) - FYIJ
               V(J4)   = V(J4) - FZIJ

!      CALCULATE CONTRIBUTION TO GORQUES

               DO K = 1, 3

                  GI(K) = VA*NR(K) + VG*X(J6-3+K)
                  GJ(K) = VB*NR(K) + VG*X(J5-3+K)

               ENDDO

               GX(J1) = GX(J1) + GI(1)
               GY(J1) = GY(J1) + GI(2)
               GZ(J1) = GZ(J1) + GI(3)
 
               GX(J2) = GX(J2) + GJ(1) 
               GY(J2) = GY(J2) + GJ(2)
               GZ(J2) = GZ(J2) + GJ(3)

            ENDDO 

         ENDDO

!     CALCULATE THE GAY-BERNE POTENTIAL 

         EGB = 4.D0 * EGB

         DO I = 1, REALNATOMS

            J      = OFFSET + 3 * I
            X(J-2) = THETA(I)
            X(J-1) = PHI(I)
            X(J)   = 0.D0
            CTI    = DCOS(THETA(I))
            STI    = DSIN(THETA(I))
            CPI    = DCOS(PHI(I))
            SPI    = DSIN(PHI(I))

            V(J-2) =  GX(I)*CTI*CPI + GY(I)*CTI*SPI - GZ(I)*STI
            V(J-1) = -GX(I)*STI*SPI + GY(I)*STI*CPI

      ENDDO

      ENDIF

      RETURN
      END SUBROUTINE GBGTP


      SUBROUTINE AATOEULER(PX,PY,PZ,PHI,THETA,CHI)

      USE COMMONS

      IMPLICIT NONE
      
      INTEGER          :: I, J, K
      DOUBLE PRECISION :: R(3), P(3), A(3), RM(3,3),PI
      DOUBLE PRECISION, INTENT(OUT) :: PHI, THETA, CHI
      DOUBLE PRECISION, INTENT(IN) :: PX, PY, PZ
 
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



!     ----------------------------------------------------------------------------------------------

      SUBROUTINE TAKESTEPELPSD (NP)

!     THIS ROUTINE TAKES STEP FOR SINGLE-SITE ELLIPSOIDAL BODIES ENSURING NO OVERLAP

      USE COMMONS 
      USE PYMODULE, ONLY : ANGLE,ANGLE2,TWOPI

      IMPLICIT NONE

      INTEGER          :: NP, JMAX, JMAX2, REALNATOMS, OFFSET, PTINDX
      INTEGER          :: J1, J2, J3, J4, J5, J6, I, J
      LOGICAL          :: OVRLPT

      DOUBLE PRECISION :: PI, DUMMY, DUMMY2, SAVECOORDS(3*NATOMS)
      DOUBLE PRECISION :: DIST(3*NATOMS/2), XMASS, YMASS, ZMASS, DMAX, VMAX, VMAX2
      DOUBLE PRECISION :: VMIN, CMMAX, CMDIST(NATOMS/2), LOCALSTEP
      DOUBLE PRECISION :: DPRAND, RANDOM, THETA, PHI
      DOUBLE PRECISION :: AE(3,3), BE(3,3) 
      DOUBLE PRECISION :: RI(3), RJ(3), RIJ(3)
      DOUBLE PRECISION :: P1(3), P2(3), ABSRIJ, RCUT, XMIN, FMIN, ECFVAL
! SF344 ADDITIONS
      LOGICAL          :: DISSOCIATED(NATOMS/2),DISSOC
      INTEGER          :: CLOSESTATOMINDEX, K1, K2, K
      DOUBLE PRECISION :: MINDISTANCE, X1, X2, Y1, Y2, Z1, Z2, RCUTARRAY(NATOMS/2), RVEC(NATOMS/2,NATOMS/2)
!      DOUBLE PRECISION, ALLOCATABLE :: RVEC(:)

!      WRITE(*,*) 'NATOMS IN TAKESTEPELPSD ', NATOMS
!      IF(ALLOCATED(RVEC)) DEALLOCATE(RVEC)

      PI         = 4.D0*ATAN(1.0D0)
      IF (GBT .OR. GBDT) THEN

         RCUT       = GBEPSNOT*MAX(GBKAPPA,1.D0)

      ELSE IF (PYGPERIODICT.OR.PYBINARYT) THEN
        DO J1=1,NATOMS/2
         RCUTARRAY(J1) = 2.0D0 * MAXVAL(PYA1BIN(J1,:)) 
        END DO
         RCUT = MAXVAL(RCUTARRAY)
      ELSE IF (LJCAPSIDT) THEN
        PYA1BIN(:,:)=0.5D0
        DO J1=1,NATOMS/2
         RCUTARRAY(J1) = 2.0D0
        END DO
         RCUT = 2.0D0
      ENDIF

      REALNATOMS = NATOMS/2
      OFFSET     = 3 * REALNATOMS

                  SAVECOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS,NP)

      DO J1 = 1,REALNATOMS

         J2     = 3*J1
        IF (PARAMONOVPBCX) THEN
!ENSURE X COMPONENT OF PARTICLE 1 VECTOR IS WITHIN BOXLX/2 OF ZERO. IF IT ISN'T THEN SUBTRACT INTEGER NUMBER OF BOXLX'S SUCH THAT IT IS.
                COORDS(J2-2,NP)=COORDS(J2-2,NP)-BOXLX*NINT(COORDS(J2-2,NP)/BOXLX)
        ENDIF

        IF (PARAMONOVPBCY) THEN
!ENSURE Y COMPONENT OF PARTICLE 1 VECTOR IS WITHIN BOXLY/2 OF ZERO. IF IT ISN'T THEN SUBTRACT INTEGER NUMBER OF BOXLY'S SUCH THAT IT IS.
                COORDS(J2-1,NP)=COORDS(J2-1,NP)-BOXLY*NINT(COORDS(J2-1,NP)/BOXLY)
        END IF

        IF (PARAMONOVPBCZ) THEN
!ENSURE Z COMPONENT OF PARTICLE 1 VECTOR IS WITHIN BOXLZ/2 OF ZERO. IF IT ISN'T THEN SUBTRACT INTEGER NUMBER OF BOXLZ'S SUCH THAT IT IS.
                COORDS(J2,NP)=COORDS(J2,NP)-BOXLZ*NINT(COORDS(J2,NP)/BOXLZ)
        ENDIF

         DUMMY2 = COORDS(J2-2,NP)**2 + COORDS(J2-1,NP)**2 + COORDS(J2,NP)**2

         IF (DUMMY2 .GT. RADIUS) THEN
! BRING BACK THE MOLECULE WITHIN THE RADIUS
                COORDS(J2-2,NP)=COORDS(J2-2,NP)-SQRT(RADIUS)*NINT(COORDS(J2-2,NP)/SQRT(RADIUS))
                COORDS(J2-1,NP)=COORDS(J2-1,NP)-SQRT(RADIUS)*NINT(COORDS(J2-1,NP)/SQRT(RADIUS))
                COORDS(J2,NP)=COORDS(J2,NP)-SQRT(RADIUS)*NINT(COORDS(J2,NP)/SQRT(RADIUS))
    WRITE(MYUNIT,'(A,2F20.10)') 'INITIAL COORDINATE OUTSIDE CONTAINER -- BRINGING MOLECULE BACK WITHIN THE CONTAINER RADIUS'

!            WRITE(*,'(A,I5,5F20.10)') 'J1,RAD,R**2,X,Y,Z:', J1, RADIUS, DUMMY2, COORDS(J2-2,NP), &
!                                       COORDS(J2-1,NP), COORDS(J2,NP)
!            PRINT*, 'INITIAL COORDINATE OUTSIDE CONTAINER -- INCREASE CONTAINER RADIUS'
!            STOP
         END IF

      END DO

      DO J1 = 1,3*NATOMS
         COORDSO(J1,NP) = COORDS(J1,NP)
      END DO

      DO J1 = 1,NATOMS/2
         VATO(J1,NP) = VAT(J1,NP)
      END DO
!     FIND THE CENTRE OF MASS

      XMASS = 0.0D0
      YMASS = 0.0D0
      ZMASS = 0.0D0

      DO J1 = 1,NATOMS/2

         XMASS = XMASS + COORDS(3*(J1-1)+1,NP)
         YMASS = YMASS + COORDS(3*(J1-1)+2,NP)
         ZMASS = ZMASS + COORDS(3*(J1-1)+3,NP)

      ENDDO

      XMASS = XMASS/(REALNATOMS)
      YMASS = YMASS/(REALNATOMS)
      ZMASS = ZMASS/(REALNATOMS)

!     FIND THE MOST WEAKLY BOUND ATOM, JMAX, THE SECOND MOST WEAKLY BOUND ATOM, JMAX2,
!     AND THE PAIR ENERGY OF THE MOST TIGHTLY BOUND ATOM, VMIN. AN ANGULAR STEP IS
!     TAKEN FOR JMAX IF ITS PAIR ENERGY IS > ASTEP*VMIN PUTTING THE ATOM AT A RADIUS OF
!     DMAX (OR CMMAX FROM CM OF THE CLUSTER).
      JMAX  = 1
      DMAX  =  1.0D0
      VMAX  = -1.0D3
      VMAX2 = -1.0D3
      VMIN  =  0.0D0
      CMMAX =  1.0D0

      DO J1 = 1, REALNATOMS

         J2 = 3*J1
         DIST(J1)   = DSQRT( COORDS(J2-2,NP)**2 + COORDS(J2-1,NP)**2 + COORDS(J2,NP)**2)
         CMDIST(J1) = DSQRT((COORDS(J2-2,NP)-XMASS)**2+(COORDS(J2-1,NP)-YMASS)**2+(COORDS(J2,NP)-ZMASS)**2)
         IF (CMDIST(J1) .GT. CMMAX) CMMAX = CMDIST(J1)
         IF (DIST(J1) .GT. DMAX) DMAX = DIST(J1)
         IF (VAT(J1,NP) .GT. VMAX) THEN
            VMAX = VAT(J1,NP)
            JMAX = J1
         ELSE IF ((VAT(J1,NP).LT. VMAX) .AND. (VAT(J1,NP) .GT. VMAX2)) THEN
              VMAX2 = VAT(J1,NP)
              JMAX2 = J1
         ENDIF
         IF (VAT(J1,NP) .LT. VMIN) VMIN = VAT(J1,NP)

      ENDDO

      IF (VAT(JMAX,NP) > (ASTEP(NP)*VMIN) .AND. (.NOT.NORESET)) THEN

         J2 = 3*JMAX
         THETA           = DPRAND()*PI
         PHI             = DPRAND()*PI*2.0D0
         COORDS(J2-2,NP) = XMASS + (CMMAX+1.0D0)*DSIN(THETA)*DCOS(PHI)
         COORDS(J2-1,NP) = YMASS + (CMMAX+1.0D0)*DSIN(THETA)*DSIN(PHI)
         COORDS(J2,NP)   = ZMASS + (CMMAX+1.0D0)*DCOS(THETA)
         DUMMY           = COORDS(J2-2,NP)**2 + COORDS(J2-1,NP)**2 + COORDS(J2,NP)**2

         IF (DUMMY > RADIUS) THEN
!          RADIUS=1.0
            DUMMY           = DSQRT(RADIUS*0.99D0/DUMMY)
            COORDS(J2-2,NP) = COORDS(J2-2,NP)*DUMMY
            COORDS(J2-1,NP) = COORDS(J2-1,NP)*DUMMY
            COORDS(J2,NP)   = COORDS(J2,NP)*DUMMY
         END IF

      ENDIF

      DO J1 = 1, REALNATOMS

         J3 = 3*J1
         J5 = OFFSET + J3

!     CHECK FOR OVERLAP

         OVRLPT = .TRUE.
!        WRITE(*,*) 'IN HERE'



95       DO WHILE (OVRLPT)
            LOCALSTEP = 0.0D0
            IF (TMOVE(NP)) LOCALSTEP = STEP(NP)

!            IF(FROZEN(J1))  LOCALSTEP = 0.0D0

            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
            COORDS(J3-2,NP) = COORDS(J3-2,NP) + LOCALSTEP*RANDOM
            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
            COORDS(J3-1,NP) = COORDS(J3-1,NP) + LOCALSTEP*RANDOM
            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
            COORDS(J3,NP)   = COORDS(J3,NP) + LOCALSTEP*RANDOM

            LOCALSTEP = 0.0D0
            IF (OMOVE(NP)) LOCALSTEP = OSTEP(NP) 

            IF(PYBINARYT) LOCALSTEP = 0.0D0

            IF(FROZEN(J1)) LOCALSTEP = 0.0D0

            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
            COORDS(J5-2,NP) = COORDS(J5-2,NP) + LOCALSTEP*RANDOM
            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
            COORDS(J5-1,NP) = COORDS(J5-1,NP) + LOCALSTEP*RANDOM
            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
            COORDS(J5,NP)   = COORDS(J5,NP) + LOCALSTEP*RANDOM
!   
            OVRLPT = .FALSE.

            RI(:)  = COORDS(J3-2:J3,NP)
            P1(:)  = COORDS(J5-2:J5,NP)

!     ROTATION MATRIX

            CALL ELPSDRTMT (P1, PYA1BIN(J1,:), AE)

            DO J2 = 1, REALNATOMS

               IF (J2 == J1) CYCLE
 
               J4    = 3*J2
               J6    = OFFSET + J4
               RJ(:) = COORDS(J4-2:J4,NP)
               P2(:) = COORDS(J6-2:J6,NP)

!            IF(FROZEN(J2))  LOCALSTEP = 0.0D0

!            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
!            COORDS(J4-2,NP) = COORDS(J4-2,NP) + LOCALSTEP*RANDOM
!            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
!            COORDS(J4-1,NP) = COORDS(J4-1,NP) + LOCALSTEP*RANDOM
!            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
!            COORDS(J4,NP)   = COORDS(J4,NP) + LOCALSTEP*RANDOM

!            LOCALSTEP = 0.0D0
!            IF (OMOVE(NP)) LOCALSTEP = OSTEP(NP) 
!            IF(PYBINARYT) LOCALSTEP = 0.0D0

!            IF(FROZEN(J2)) LOCALSTEP = 0.0D0

!            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
!            COORDS(J6-2,NP) = COORDS(J6-2,NP) + LOCALSTEP*RANDOM
!            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
!            COORDS(J6-1,NP) = COORDS(J6-1,NP) + LOCALSTEP*RANDOM
!            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
!            COORDS(J6,NP)   = COORDS(J6,NP) + LOCALSTEP*RANDOM

!     ROTATION MATRIX

               CALL ELPSDRTMT (P2, PYA1BIN(J2,:), BE)

               RIJ    = RI - RJ
               ABSRIJ = DSQRT(DOT_PRODUCT(RIJ,RIJ))
                        
               IF (ABSRIJ < RCUT) THEN 

!     DETERMINE ELLIPTIC CONTACT FUNCTION

                  CALL BRENTMIN (0.D0, 0.51D0, 1.D0, AE, BE, RIJ, XMIN, FMIN)

                  ECFVAL = - FMIN
! ALLOW FOR A SLIGHT OVERLAP 
                  IF (ECFVAL < PYOVERLAPTHRESH) THEN
!                        WRITE(*,*) 'ATOMS OVERLAPPING', J1, J2, ECFVAL, ABSRIJ
                     OVRLPT = .TRUE.
                     GO TO 95
                  ENDIF

               ENDIF

            ENDDO  ! END LOOP OVER J2
         ENDDO  ! END WHILE
      ENDDO  ! END LOOP OVER J1   

                     IF (FREEZE) THEN
                        DO J2=1,NATOMS
                           IF (FROZEN(J2)) THEN
                              COORDS(3*(J2-1)+1:3*(J2-1)+3,NP)=SAVECOORDS(3*(J2-1)+1:3*(J2-1)+3)
                           ENDIF
                        ENDDO
                     ENDIF


! NOW DETERMINE WHETHER ANY PARTICLE HAS DISSOCIATED
 IF(PYGPERIODICT) THEN
  DO
  DISSOC=.FALSE.
   DO J1=1,NATOMS/2-1
    J2=3*J1
    X1=COORDS(J2-2,NP)
    Y1=COORDS(J2-1,NP)
    Z1=COORDS(J2,NP)
   
    DISSOCIATED(J1)=.TRUE.
    DO K1=1,NATOMS/2
      K2=3*K1
      X2=COORDS(K2-2,NP)
      Y2=COORDS(K2-1,NP)
      Z2=COORDS(K2,NP)
      RVEC(J1,K1)=SQRT((X2-X1)**2+(Y2-Y1)**2+(Z2-Z1)**2)
!      WRITE(*,*) J1,K1,RVEC(J1,K1), 4*RCUT
      IF(RVEC(J1,K1)<4*RCUT) DISSOCIATED(J1)=.FALSE.
!      WRITE(*,*) 'IN HERE, CUTOFF=',PCUTOFF,J1
    END DO
    IF(DISSOCIATED(J1)) THEN
        ! DETERMINE THE CLOSEST NEIGHBOR
        MINDISTANCE = HUGE(1.0D0)
        DO K1=1,NATOMS/2
         IF(K1==J1) CYCLE
         IF(RVEC(J1,K1)<MINDISTANCE) THEN
           MINDISTANCE=RVEC(J1,K1)
           CLOSESTATOMINDEX=K1
         END IF
        END DO
        WRITE(MYUNIT,'(A,I6,A,F10.6,I6)') 'ATOM ', J1, ' IS DISSOCIATED, DISTANCE TO THE CLOSEST ATOM: ', &
        & MINDISTANCE, CLOSESTATOMINDEX
        DISSOC=.TRUE.
        ! MOVE THE DISSOCIATED ATOM CLOSER TO THE CLOSEST ATOM
        COORDS(J2-2,NP)=0.5D0*(COORDS(3*CLOSESTATOMINDEX-2,NP)+COORDS(J2-2,NP))
        COORDS(J2-1,NP)=0.5D0*(COORDS(3*CLOSESTATOMINDEX-1,NP)+COORDS(J2-1,NP))
        COORDS(J2  ,NP)=0.5D0*(COORDS(3*CLOSESTATOMINDEX  ,NP)+COORDS(J2  ,NP))
    END IF
   END DO
   IF(.NOT.DISSOC) EXIT
  END DO
        DO J1=NATOMS/2+1,NATOMS
            J2=3*J1
            ANGLE2=DOT_PRODUCT(COORDS(J2-2:J2,NP),COORDS(J2-2:J2,NP))

            IF(ANGLE2>TWOPI**2) THEN
                ANGLE2=SQRT(ANGLE2)
                ANGLE = ANGLE2 - INT(ANGLE2/TWOPI)*TWOPI
                COORDS(J2-2:J2,NP)=COORDS(J2-2:J2,NP)/ANGLE2 * ANGLE
            END IF
        END DO

RETURN
!     CHECK FOR OVERLAP AGAIN (NOT WORKING!)
      DO J1 = 1, REALNATOMS

         J3 = 3*J1
         J5 = OFFSET + J3

         OVRLPT = .TRUE.

195       DO WHILE (OVRLPT)

            LOCALSTEP = 0.0D0
            IF (TMOVE(NP)) LOCALSTEP = STEP(NP)

            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
            COORDS(J3-2,NP) = COORDS(J3-2,NP) + LOCALSTEP*RANDOM
            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
            COORDS(J3-1,NP) = COORDS(J3-1,NP) + LOCALSTEP*RANDOM
            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
            COORDS(J3,NP)   = COORDS(J3,NP) + LOCALSTEP*RANDOM

            LOCALSTEP = 0.0D0
            IF (OMOVE(NP)) LOCALSTEP = OSTEP(NP)

            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
            COORDS(J5-2,NP) = COORDS(J5-2,NP) + LOCALSTEP*RANDOM
            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
            COORDS(J5-1,NP) = COORDS(J5-1,NP) + LOCALSTEP*RANDOM
            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
            COORDS(J5,NP)   = COORDS(J5,NP) + LOCALSTEP*RANDOM
          
            OVRLPT = .FALSE.

            RI(:)  = COORDS(J3-2:J3,NP)
            P1(:)  = COORDS(J5-2:J5,NP)

!     ROTATION MATRIX

            CALL ELPSDRTMT (P1, ESA, AE)

            DO J2 = 1, REALNATOMS

               IF (J2 == J1) CYCLE
 
               J4    = 3*J2
               J6    = OFFSET + J4
               RJ(:) = COORDS(J4-2:J4,NP)
               P2(:) = COORDS(J6-2:J6,NP)

!     ROTATION MATRIX

               CALL ELPSDRTMT (P2, ESA, BE)

               RIJ    = RI - RJ
               ABSRIJ = DSQRT(DOT_PRODUCT(RIJ,RIJ))
                        
               IF (ABSRIJ < RCUT) THEN 

!     DETERMINE ELLIPTIC CONTACT FUNCTION

                  CALL BRENTMIN (0.D0, 0.51D0, 1.D0, AE, BE, RIJ, XMIN, FMIN)

                  ECFVAL = - FMIN
 
                  IF (ECFVAL < 1.D0) THEN
                     OVRLPT = .TRUE.
                     GO TO 195
                  ENDIF

               ENDIF

            ENDDO  ! END LOOP OVER J2

         ENDDO  ! END WHILE
 
      ENDDO  ! END LOOP OVER J1   

 END IF

      END SUBROUTINE TAKESTEPELPSD

SUBROUTINE TAKESTEPSWAPMOVES(NP)
USE COMMONS, ONLY : NATOMS, COORDS, PYBINARYTYPE1,MYUNIT,SWAPMOVEST,PYSWAP

IMPLICIT NONE

INTEGER :: NP, RANDOM1, RANDOM2, J1, J2
DOUBLE PRECISION :: COORDSSTORE(3,4),DPRAND

! GET TWO RANDOM INTEGER NUMBERS, ONE FOR ATOM INDEX TYPE 1, ONE FOR ATOM INDEX TYPE 2

DO J2=1,PYSWAP(3)
        DO
          RANDOM1=INT(DPRAND()*1000)
!           RANDOM1=J2
          IF(RANDOM1<=PYBINARYTYPE1.AND.RANDOM1>=1) EXIT
        END DO

        DO
          RANDOM2=INT(DPRAND()*1000)
          IF(RANDOM2>PYBINARYTYPE1.AND.RANDOM2<=NATOMS/2) EXIT
        END DO
        !RANDOM1=PYSWAP(1)
        !RANDOM2=PYSWAP(2)
        COORDSSTORE(:,:)=0.0D0
        ! SELECT A BODY FROM THE TWO TYPES
        DO J1=1,3
         COORDSSTORE(J1,1)=COORDS(3*(RANDOM1-1)+J1,NP)
         COORDSSTORE(J1,2)=COORDS(3*(RANDOM2-1)+J1,NP)
         COORDSSTORE(J1,3)=COORDS(3*NATOMS/2+3*(RANDOM1-1)+J1,NP)
         COORDSSTORE(J1,4)=COORDS(3*NATOMS/2+3*(RANDOM2-1)+J1,NP)
        END DO

        ! SWAP COORDINATES
        WRITE(MYUNIT,*) 'SWAPPING ATOMS ', RANDOM1, RANDOM2
        !WRITE(MYUNIT,*) COORDSSTORE(:,1)
        !WRITE(MYUNIT,*) COORDSSTORE(:,3)
        !WRITE(MYUNIT,*) COORDSSTORE(:,2)
        !WRITE(MYUNIT,*) COORDSSTORE(:,4)

        DO J1=1,3
         COORDS(3*(RANDOM1-1)+J1,NP)=COORDSSTORE(J1,2)
         COORDS(3*(RANDOM2-1)+J1,NP)=COORDSSTORE(J1,1)
         COORDS(3*NATOMS/2+3*(RANDOM1-1)+J1,NP)=COORDSSTORE(J1,4)
         COORDS(3*NATOMS/2+3*(RANDOM2-1)+J1,NP)=COORDSSTORE(J1,3)
        END DO
        
        IF(PYSWAP(2)<NATOMS/2) THEN 
                PYSWAP(2)=PYSWAP(2)+1
        ELSE
                PYSWAP(1)=PYSWAP(1)+1
                PYSWAP(2)=PYBINARYTYPE1+1
        END IF
        IF(PYSWAP(1)>PYBINARYTYPE1) THEN 
              WRITE(MYUNIT,*) 'ALL ATOMS SWAPPED, RESTARTING'
!      SWAPMOVEST=.FALSE.  
              PYSWAP(1)=1
              PYSWAP(2)=PYBINARYTYPE1+1
              EXIT
        END IF
END DO
END SUBROUTINE TAKESTEPSWAPMOVES

SUBROUTINE PYPES(XORIG)
USE COMMONS,ONLY : NATOMS
IMPLICIT NONE

INTEGER :: J0,J1,J2,J3,NSTEPS
DOUBLE PRECISION :: INTERVAL, XORIG(3*NATOMS),XNEW(3*NATOMS),GRAD(3*NATOMS)
DOUBLE PRECISION, ALLOCATABLE   :: PES(:,:,:),X1(:),Y1(:),Z1(:)

INTERVAL=0.2D0
NSTEPS=60

XNEW(:)=XORIG(:)

ALLOCATE(PES(NSTEPS+1,NSTEPS+1,NSTEPS+1),X1(NSTEPS+1),Y1(NSTEPS+1),Z1(NSTEPS+1))
OPEN(UNIT=955,FILE='PY_PES',STATUS='UNKNOWN')
DO J0=0,NSTEPS
 XNEW(2)=XORIG(2)+J0*INTERVAL
 Y1(J0+1)=XNEW(2)
 DO J1=0,NSTEPS
   WRITE(*,*) 'WORKING... ', J0,J1
   XNEW(1)=XORIG(1)+J1*INTERVAL
   X1(J1+1)=XNEW(1)
   DO J2=0,NSTEPS
      XNEW(3)=XORIG(3)+J2*INTERVAL
      Z1(J2+1)=XNEW(3)
      CALL PYGPERIODIC(XNEW,GRAD,PES(J0+1,J1+1,J2+1),.TRUE.)
      IF(PES(J0+1,J1+1,J2+1)>1.0D9) THEN
        PES=1.0D9
      ELSE IF(PES(J0+1,J1+1,J2+1)<-1.0D3) THEN
        PES=-1.0D9
      END IF
         WRITE(955,'(4F20.5)') X1(J1+1),Y1(J0+1),Z1(J2+1),PES(J0+1,J1+1,J2+1)
   END DO
 END DO
END DO

DO J0=1,NSTEPS+1
 DO J1=1,NSTEPS+1
   DO J2=1,NSTEPS+1
!        WRITE(955,'(4F20.5)') X1(J1),Y1(J0),Z1(J2),PES(J0,J1,J2)
   END DO
 END DO
END DO
CLOSE(955)
END SUBROUTINE PYPES

SUBROUTINE INITIALISELJCAPSIDMODEL

USE COMMONS, ONLY: BOXLX,BOXLY,BOXLZ,PARAMONOVPBCX,PARAMONOVPBCY,PARAMONOVPBCZ,PCUTOFF,PARAMONOVCUTOFF,&
                &       NATOMS,PYA1BIN,PYA2BIN,PYSIGNOT,PYEPSNOT,RADIFT,LJSITE,BLJSITE,PEPSILON1,&
                &       PSCALEFAC1,PSCALEFAC2,MAXINTERACTIONS,PYBINARYT,PYBINARYTYPE1,MYUNIT,VT
USE LJCAPSIDMODULE 

IMPLICIT NONE
   WRITE(MYUNIT,*) 'INITIALISING VARIABLES FOR LJ CAPSID MODEL',NATOMS 
! ALLOCATE ARRAYS

    ALLOCATE(RMIVEC(NATOMS/2,3,3),DPI1RMVEC(NATOMS/2,3,3), DPI2RMVEC(NATOMS/2,3,3), DPI3RMVEC(NATOMS/2,3,3))
    ALLOCATE(PSCALEFAC1VEC(NATOMS/2),PSCALEFAC2VEC(NATOMS/2))
    ALLOCATE(AEZR1(NATOMS/2,3,3), AEZR2(NATOMS/2,3,3))
    ALLOCATE(EPSILON1(4,NATOMS/2,NATOMS/2))
      I3(:,:)    = 0.D0
      AEZR1(:,:,:) = 0.D0
      AEZR2(:,:,:) = 0.D0



     PI=ATAN(1.0D0)*4.0
     TWOPI=2.0D0*PI

       EPSILON1(1:3,:,:)=PEPSILON1(1)   ! WE ARE GOING TO USE EPSILON1 FOR THE EXTRA LJ SITES
       EPSILON1(4,:,:)=PYEPSNOT
       SIGMA1=1.0D0 ! PSIGMA1 IS NONEXISTENT FROM NOW ON
       PSCALEFAC1VEC(:)=PSCALEFAC1(1)
       PSCALEFAC2VEC(:)=PSCALEFAC2(1)

      REALNATOMS = NATOMS/2
      OFFSET     = 3*REALNATOMS
      
      DO K1 = 1, 3
        I3(K1,K1) = 1.0D0
      ENDDO
END SUBROUTINE INITIALISELJCAPSIDMODEL

SUBROUTINE LJCAPSIDMODEL (X, G, ENERGY, GTEST)
USE COMMONS, ONLY:  NATOMS,PYSIGNOT,PYEPSNOT,LJSITE,BLJSITE,PEPSILON1,&
                &       PSCALEFAC1,PSCALEFAC2,MAXINTERACTIONS,PYBINARYT,PYBINARYTYPE1,MYUNIT,VT
USE LJCAPSIDMODULE

IMPLICIT NONE
DOUBLE PRECISION :: DUMMY, VDUMMY, DVDUMMY(12), DUMMY1, DUMMY2, P(3), TERM2(MAXINTERACTIONS), TERM3(MAXINTERACTIONS), &
                  & XLJ(MAXINTERACTIONS,2,3), RLJVEC(MAXINTERACTIONS,3), RLJUNITVEC(MAXINTERACTIONS,3), VECSBF(3),&
                  & DRLJ(MAXINTERACTIONS,12), RLJ(MAXINTERACTIONS), RLJ2(MAXINTERACTIONS), &
                  & LLJ(12,MAXINTERACTIONS), DLLJ1(MAXINTERACTIONS,12),ATTR(MAXINTERACTIONS), &
                  & TIJ(3), TJI(3), FIJ(3), ENERGY, X(3*NATOMS), G(3*NATOMS), SIGMA1VEC(MAXINTERACTIONS), &
                  & RMI(3,3), RMJ(3,3), RI(3), RJ(3), &
                  & DPI1RM(3,3), DPI2RM(3,3), DPI3RM(3,3), DPJ1RM(3,3), DPJ2RM(3,3), DPJ3RM(3,3)
INTEGER          :: K
!INTEGER          :: J1, J3, J5, K, K1, REALNATOMS
LOGICAL          :: GTEST


VT(1:NATOMS/2)=0.0D0
 
TERM2(:)=1.0D0
TERM3(:)=0.0D0
ATTR(1:3)=0.0D0
ATTR(4)=1.0D0
SIGMA1VEC(1:3)=1.0
SIGMA1VEC(4)=1.0
!PI=ATAN(1.0D0)*4.0
!TWOPI=2.0D0*PI

!      REALNATOMS = NATOMS/2
!      OFFSET     = 3*REALNATOMS

!      DO K1 = 1, 3
!        I3(K1,K1) = 1.0D0
!      ENDDO

      ENERGY = 0.D0
      G(:)   = 0.D0

        DO J1=1, REALNATOMS
            J3      = 3*J1
            J5      = OFFSET + J3
            RI      = X(J3-2:J3)
            P       = X(J5-2:J5)
            ANGLE=SQRT(DOT_PRODUCT(P,P))
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
!        WRITE(*,*) 'CALLING RMDRVT', P, RMIVEC(J1,:,:), DPI1RMVEC(J1,:,:), DPI2RMVEC(J1,:,:), DPI3RMVEC(J1,:,:), GTEST
            CALL RMDRVT(P, RMIVEC(J1,:,:), DPI1RMVEC(J1,:,:), DPI2RMVEC(J1,:,:), DPI3RMVEC(J1,:,:), GTEST)

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


           DO J2 = J1 + 1, REALNATOMS
               J4     = 3*J2
               J6     = OFFSET + J4
               RJ     = X(J4-2:J4)
               P      = X(J6-2:J6)

               RMJ(:,:)=RMIVEC(J2,:,:)
               DPJ1RM(:,:)=DPI1RMVEC(J2,:,:)
               DPJ2RM(:,:)=DPI2RMVEC(J2,:,:)
               DPJ3RM(:,:)=DPI3RMVEC(J2,:,:)

               VDUMMY=0.0D0
               DVDUMMY(:)=0.0D0
            DO K=1,MAXINTERACTIONS ! K=1 -- INTERACTION BETWEEN REPULSIVE PRIMARY 'APEX' SITES
                         ! K=2 AND K=3 -- INTERACTION BETWEEN SECONDARY AND PRIMARY 'APEX' SITES 
                         ! K=4  -- LJ SITE AT THE CENTRE (NORMAL 12-6 INTERACTION)
                        VECSBF(1)=1.0D0
                        VECSBF(2)=0.0D0
                        VECSBF(3)=0.0D0
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
                   ELSE IF(K==4) THEN
                        DUMMY1=0.0D0
                        DUMMY2=0.0D0
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

              ! INTERACTION BETWEEN THE LJ SITES:
                        
                        LLJ(1,K)=SIGMA1VEC(K)*DUMMY !/RLJ(K)
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
                                DLLJ1(K,J) =-SIGMA1VEC(K)*DUMMY*DUMMY*DRLJ(K,J)
                            END DO
                 VDUMMY=VDUMMY+4.0D0*EPSILON1(K,J1,J2)*TERM2(K)*(LLJ(12,K) - ATTR(K)*LLJ(6,K))
                                     ! EXTRA LJ SITES, REPULSIVE, AND LJ SITE IN THE MIDDLE (ATTR(4)=1, ATTR(1:3)=0)
!                 IF(K==3) WRITE(*,*) 'REPULSIVE ENERGIES ', 4.0D0*EPSILON1(K,J1,J2)*LLJ(12,K)*TERM2(K)
             DO J=1,3
               DVDUMMY(J) = DVDUMMY(J) + 4.0D0*EPSILON1(K,J1,J2)*(12.0D0*LLJ(11,K)*DLLJ1(K,J))*TERM2(K) + &
                          & 4.0D0*EPSILON1(K,J1,J2)*LLJ(12,K)*TERM3(K)*DRLJ(K,J)! EXTRA LJ SITE DERIVATIVES, CURRENTLY ONLY REPULSIVE
               DVDUMMY(J) = DVDUMMY(J) - ATTR(K)*(4.0D0*EPSILON1(K,J1,J2)*(6.0D0*LLJ(5,K)*DLLJ1(K,J))*TERM2(K) + &
                          & 4.0D0*EPSILON1(K,J1,J2)*LLJ(6,K)*TERM3(K)*DRLJ(K,J))
               FIJ(J) = DVDUMMY(J)
             END DO
             DO J=7,12
               DVDUMMY(J) = DVDUMMY(J) + 4.0D0*EPSILON1(K,J1,J2)*(12.0D0*LLJ(11,K)*DLLJ1(K,J))*TERM2(K) + &
                          & 4.0D0*EPSILON1(K,J1,J2)*LLJ(12,K)*TERM3(K)*DRLJ(K,J)! EXTRA LJ SITE DERIVATIVES, CURRENTLY ONLY REPULSIVE
               DVDUMMY(J) = DVDUMMY(J) - ATTR(K)*(4.0D0*EPSILON1(K,J1,J2)*(6.0D0*LLJ(5,K)*DLLJ1(K,J))*TERM2(K) + &
                          & 4.0D0*EPSILON1(K,J1,J2)*LLJ(6,K)*TERM3(K)*DRLJ(K,J))
             END DO
            END DO !K=1,4

               TIJ = DVDUMMY(7:9)
               TJI = DVDUMMY(10:12)

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

               ENERGY = ENERGY + VDUMMY
               VT(J1) = VT(J1) + VDUMMY
               VT(J2) = VT(J2) + VDUMMY        ! PAIR POTENTIALS

               G(J3-2:J3) = G(J3-2:J3) - FIJ
               G(J4-2:J4) = G(J4-2:J4) + FIJ
               G(J5-2:J5) = G(J5-2:J5) + TIJ
               G(J6-2:J6) = G(J6-2:J6) + TJI

        END DO
       END DO

END SUBROUTINE LJCAPSIDMODEL

