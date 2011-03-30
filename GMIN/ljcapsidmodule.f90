
MODULE LJCAPSIDMODULE

 INTEGER          :: I, J, J1, J2, J3, J4, J5, J6, K1, K2, OFFSET, REALNATOMS
 DOUBLE PRECISION, ALLOCATABLE :: RMIvec(:,:,:), DPI1RMvec(:,:,:), DPI2RMvec(:,:,:), DPI3RMvec(:,:,:)
 DOUBLE PRECISION, ALLOCATABLE :: PSCALEFAC1vec(:),PSCALEFAC2vec(:),epsilon1(:,:,:),AEZR1(:,:,:), AEZR2(:,:,:)

 DOUBLE PRECISION :: angle,angle2,pi,twopi,sigma1,cut,ron,ron2,range2inv3

 DOUBLE PRECISION ::I3(3,3) 

END MODULE LJCAPSIDMODULE

SUBROUTINE INITIALISEPYGPERIODIC

use commons, only: BOXLX,BOXLY,BOXLZ,PARAMONOVPBCX,PARAMONOVPBCY,PARAMONOVPBCZ,PCUTOFF,PARAMONOVCUTOFF,&
                &       natoms,pya1bin,pya2bin,pysignot,pyepsnot,radift,LJSITE,BLJSITE,PEPSILON1,&
                &       PSCALEFAC1,PSCALEFAC2,MAXINTERACTIONS,PYBINARYT,PYBINARYTYPE1,MYUNIT,VT, &
                &       PEPSILONATTR, PSIGMAATTR, LJSITEATTR, LJSITECOORDST, LJSITECOORDS

use pymodule

implicit none
   WRITE(MYUNIT,*) 'initialising variables for PY',NATOMS 
! allocate arrays

    ALLOCATE(RMIvec(natoms/2,3,3),DPI1RMvec(natoms/2,3,3), DPI2RMvec(natoms/2,3,3), DPI3RMvec(natoms/2,3,3))
    ALLOCATE(PSCALEFAC1vec(natoms/2),PSCALEFAC2vec(natoms/2),epsilon1(4,natoms/2,natoms/2))
    ALLOCATE(AEZR1(NATOMS/2,3,3), AEZR2(NATOMS/2,3,3))
    IF(.NOT.ALLOCATED(VT)) ALLOCATE(VT(NATOMS/2))
    IF(.NOT.ALLOCATED(PYA1bin)) ALLOCATE(PYA1bin(NATOMS/2,3))
    IF(.NOT.ALLOCATED(PYA2bin)) ALLOCATE(PYA2bin(NATOMS/2,3))

    
          vecsbf(1)=1.0D0
          vecsbf(2)=0.0D0
          vecsbf(3)=0.0D0

    IF(LJSITECOORDST) THEN
        vecsbf(:)=LJSITECOORDS(:)
        WRITE(MYUNIT,'(A,3F8.3)') 'repulsive LJ site coordinates will be ', LJSITECOORDS(:)
    END IF
      I3(:,:)    = 0.D0
      AEZR1(:,:,:) = 0.D0
      AEZR2(:,:,:) = 0.D0

     pi=ATAN(1.0D0)*4.0
     twopi=2.0D0*pi

     IF(LJSITE) THEN
       epsilon1(:,:,:)=PEPSILON1(1)   ! we are going to use epsilon1 for the extra LJ sites
       sigma1(:)=1.0D0 ! PSIGMA1 is nonexistent from now on, except for the attractive secondary apex sites
       PSCALEFAC1vec(:)=PSCALEFAC1(1)
       PSCALEFAC2vec(:)=PSCALEFAC2(1)
       IF(LJSITEATTR) THEN !attractive secondary apex site is turned on
        attr(:)=1.0D0
        epsilon1(1,:,:)=PEPSILONATTR(1)
        epsilon1(2,:,:)=0.0D0!SQRT(PEPSILONATTR(1)*PEPSILONATTR(2))
        epsilon1(3,:,:)=epsilon1(2,:,:)
        epsilon1(4,:,:)=PEPSILONATTR(2)
        
        sigma1(1)=PSIGMAATTR(1)
        sigma1(2)=0.5D0*(PSIGMAATTR(1)+PSIGMAATTR(2))
        sigma1(3)=sigma1(2)
        sigma1(4)=PSIGMAATTR(2)
       END IF
     ELSE
       epsilon1(:,:,:)=0.0D0
       sigma1(:)=0.0D0
     END IF

! sanity checks
       IF(PYBINARYT.AND.LJSITE.AND..NOT.BLJSITE) THEN
        WRITE(MYUNIT,*) 'ERROR --- for binary PY systems with extra LJ sites '// &
                        & 'you have to specify the parameters for both types separately (>3 arguments after EXTRALJSITE)! '
        STOP
       END IF
       IF(BLJSITE.AND..NOT.PYBINARYT) THEN
        WRITE(MYUNIT,*) 'ERROR --- binary LJ sites specified, but no binary PY particles. '// &
                        & 'EXTRALJSITE should have only 3 arguments!'
        STOP
       END IF


       IF(PYBINARYT.AND.BLJSITE) THEN
        DO J1=1,NATOMS/2
          IF(J1<=PYBINARYTYPE1) THEN
                PSCALEFAC1vec(J1)=PSCALEFAC1(1)
                PSCALEFAC2vec(J1)=PSCALEFAC2(1)
          ELSE
                PSCALEFAC1vec(J1)=PSCALEFAC1(2)
                PSCALEFAC2vec(J1)=PSCALEFAC2(2)
          END IF
        END DO

        DO J1=1,NATOMS/2-1
          DO J2=J1+1,NATOMS/2
           IF(J1<=PYBINARYTYPE1.AND.J2<=PYBINARYTYPE1) THEN
                epsilon1(:,J1,J2)=PEPSILON1(1)
           ELSE IF(J1<=PYBINARYTYPE1.AND.J2>PYBINARYTYPE1) THEN
                epsilon1(:,J1,J2)=PEPSILON1(3)
           ELSE
                epsilon1(:,J1,J2)=PEPSILON1(2)
           END IF
          END DO
        END DO
       END IF
! cutoff stuff for the extra LJ sites
       cut = PCUTOFF*PCUTOFF ! cutoff squared
       ron = PCUTOFF*0.9D0
       ron2 = ron*ron
       range2inv3=1.0D0/(cut-ron2)**3
!     FROM INPUT PARAMETERS

      REALNATOMS = NATOMS/2
      OFFSET     = 3*REALNATOMS
      
      DO K1 = 1, 3
        I3(K1,K1) = 1.0D0
      ENDDO
      DO J1=1,REALNATOMS
       DO K1 = 1, 3
         AEZR1(J1,K1,K1) = 1.D0/(PYA1bin(J1,K1)*PYA1bin(J1,K1))
         AEZR2(J1,K1,K1) = 1.D0/(PYA2bin(J1,K1)*PYA2bin(J1,K1))
       ENDDO
      END DO

END SUBROUTINE INITIALISEPYGPERIODIC

! PY potential, dc430's implementation
! with PBC and continuous cutoff added
! plus extra LJ site
      SUBROUTINE PYGPERIODIC (X, G, ENERGY, GTEST)

       use commons, only: BOXLX,BOXLY,BOXLZ,PARAMONOVPBCX,PARAMONOVPBCY,PARAMONOVPBCZ,PCUTOFF,PARAMONOVCUTOFF,&
                &       natoms,pya1bin,pya2bin,pysignot,pyepsnot,radift,LJSITE,BLJSITE,PEPSILON1,& 
                &       PSCALEFAC1,PSCALEFAC2,MAXINTERACTIONS,PYBINARYT,PYBINARYTYPE1,MYUNIT,VT,PEPSILONATTR,PSIGMAATTR,&
                &       FROZEN

       use pymodule
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
      double precision :: D1ABEZ(3,3), D2ABEZ(3,3), D3ABEZ(3,3), D1ABE(3,3), D2ABE(3,3), D3ABE(3,3) 
      LOGICAL          :: GTEST


! sf344 additions

      DOUBLE PRECISION :: CLJ1(2),CLJ2(2),CLJ3(2),CLJ4(2),CLJ5(2),CLJ6(2),CLJ7(2),CLJ8(2),CLJ11(2),CLJ12(2), &
                          & CLJ13(2),CLJ14(2),DFDR(2,3),DFP(2,6),dr(6),dCLJ1(2,12),dVDUMMY(12),LJ1(2),DUMMY,DUMMY1,DUMMY2,&
                          & dLJ1(2,12),VDUMMY
      DOUBLE PRECISION :: term2(maxinteractions),term3(maxinteractions), &
                          & xlj(maxinteractions,2,3),rljvec(maxinteractions,3),rljunitvec(maxinteractions,3),&
                          & drlj(maxinteractions,12),rlj(maxinteractions),rlj2(maxinteractions)
      DOUBLE PRECISION :: LLJ(12,maxinteractions), dLLJ1(maxinteractions,12)
      INTEGER          :: k

     VT(1:NATOMS/2)=0.0D0
 
     term2(:)=1.0D0
     term3(:)=0.0D0
      
!       IF(PYBINARYT) THEN
!        DO J1=1,NATOMS/2-1
!          DO J2=J1+1,NATOMS/2
!           IF(J1<=PYBINARYTYPE1.AND.J2>PYBINARYTYPE1) THEN
!                epsilon1(J1,J2)=PEPSILON1/10.0D0
!           ELSE
!                epsilon1(J1,J2)=PEPSILON1
!          END DO
!        END DO
!       END IF

      DO J1=1,REALNATOMS
        J2=3*J1

        IF (PARAMONOVPBCX) THEN
!ensure x component of particle 1 vector is within BoxLx/2 of zero. If it isn't then subtract integer number of boxlx's such that it is.
                X(J2-2)=X(J2-2)-BOXLX*NINT(X(J2-2)/BOXLX)
        ENDIF

        IF (PARAMONOVPBCY) THEN
!ensure y component of particle 1 vector is within BoxLy/2 of zero. If it isn't then subtract integer number of boxly's such that it is.
                X(J2-1)=X(J2-1)-BOXLY*NINT(X(J2-1)/BOXLY)
        END IF

        IF (PARAMONOVPBCZ) THEN
!ensure z component of particle 1 vector is within BoxLz/2 of zero. If it isn't then subtract integer number of boxlz's such that it is.
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
            angle=sqrt(dot_product(P,P))
            if(angle>twopi) then
! normalise angle-axis coordinates
                X(J5-2:J5)=X(J5-2:J5)/angle
                do
                  angle=angle-twopi
                  if(angle<2*pi) exit
                end do
! multiply with new angle
                X(J5-2:J5)=X(J5-2:J5)*angle
            end if

            CALL RMDRVT(P, RMIvec(J1,:,:), DPI1RMvec(J1,:,:), DPI2RMvec(J1,:,:), DPI3RMvec(J1,:,:), GTEST)

        END DO

         DO J1 = 1, REALNATOMS
          
            J3      = 3*J1
            J5      = OFFSET + J3
            RI      = X(J3-2:J3)
            P       = X(J5-2:J5)
!     ROTATION MATRIX

!            CALL RMDRVT(P, RMI, DPI1RM, DPI2RM, DPI3RM, GTEST)
            RMI(:,:)=RMIvec(J1,:,:)
            DPI1RM(:,:)=DPI1RMvec(J1,:,:)
            DPI2RM(:,:)=DPI2RMvec(J1,:,:)
            DPI3RM(:,:)=DPI3RMvec(J1,:,:)

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
               RMJ(:,:)=RMIvec(J2,:,:)
               DPJ1RM(:,:)=DPI1RMvec(J2,:,:)
               DPJ2RM(:,:)=DPI2RMvec(J2,:,:)
               DPJ3RM(:,:)=DPI3RMvec(J2,:,:)
     
               BE1 = MATMUL(RMJ,(MATMUL(AEZR1(J2,:,:),(TRANSPOSE(RMJ)))))

               IF (RADIFT) THEN
   
                  BE2 = MATMUL(RMJ,(MATMUL(AEZR2(J2,:,:),(TRANSPOSE(RMJ)))))

               ENDIF

!     CALCULATE SEPARATION

               RIJ    = RI - RJ
               RIJSQ  = DOT_PRODUCT(RIJ,RIJ)
               ABSRIJ = DSQRT(RIJSQ)
               NR     = RIJ / ABSRIJ

                IF(PARAMONOVCUTOFF.AND.RIJSQ>cut) GOTO 124

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


! Correction terms to the potential if we require a cutoff at rc. 
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
! r/SQRT(F1(k)) = f(A,B), a parameter ONLY dependent on orientation!
                         dr(1)=1.0D0/ABSRIJ*(RI(1)-RJ(1))
                         dr(2)=1.0D0/ABSRIJ*(RI(2)-RJ(2))
                         dr(3)=1.0D0/ABSRIJ*(RI(3)-RJ(3))
                         dr(4)=-dr(1)
                         dr(5)=-dr(2)
                         dr(6)=-dr(3)
                     do k=1,2
                          CLJ1(k)=PYSIGNOT/(PCUTOFF-ABSRIJ*SRTFI(k)+PYSIGNOT)
                          CLJ2(k)=CLJ1(k)**2
                          CLJ3(k)=CLJ2(k)*CLJ1(k)
                          CLJ4(k)=CLJ2(k)**2
                          CLJ5(k)=CLJ4(k)*CLJ1(k)
                          CLJ6(k)=CLJ4(k)*CLJ2(k)
                          CLJ7(k)=CLJ6(k)*CLJ1(k)
                          CLJ8(k)=CLJ6(k)*CLJ2(k)
                          CLJ11(k)=CLJ5(k)*CLJ6(k)
                          CLJ12(k)=CLJ6(k)**2
                          CLJ13(k)=CLJ12(k)*CLJ1(k)
                          CLJ14(k)=CLJ7(k)**2

                         DUMMY=CLJ1(k)/PYSIGNOT
                         DUMMY=DUMMY**2
                         DUMMY1=SRTFI(k)
                         DUMMY2=DUMMY1**3

                         DO j=1,3
                            dCLJ1(k,j) = 1.0D0*PYSIGNOT*DUMMY*(-1.0D0*DUMMY1*dr(j)+0.5D0*ABSRIJ*DUMMY2*DFDR(k,j))
                            dCLJ1(k,j+3)= -dCLJ1(k,j)
                         END DO
                         DO j=7,12
                            dCLJ1(k,j) =-1.0D0*PYSIGNOT*DUMMY*(0.5D0*ABSRIJ*DUMMY2*DFP(k,j-6)) !derivatives wrt to orientation
                         END DO
                     end do
                   END IF
 !     CALCULATE PY POTENTIAL ENERGY
123 CONTINUE 
        VDUMMY = 0.0D0
                 
!        IF(FROZEN(J1)) THEN 
                VDUMMY=0.0D0
!        ELSE
               IF(LJSITE) THEN  
                do k=1,maxinteractions ! k=1 -- interaction between repulsive primary 'apex' sites
                         ! k=2 and k=3 -- interaction between secondary and primary 'apex' sites
                         ! k=4 -- interaction between secondary 'apex' sites (normal LJ interaction) 
                        ! trying to modify code to allow for binary systems. 
                        ! apex site heights will be defined in absolute units,
                        ! hence PYA1bin(J1,1) etc. will be removed from below

                   IF(k==1) THEN
                        DUMMY1=PSCALEFAC1vec(J1)!*PYA1bin(J1,1)
                        DUMMY2=PSCALEFAC1vec(J2)!*PYA1bin(J2,1)
                   ELSE IF(k==2) THEN
                        DUMMY1=PSCALEFAC1vec(J1)!*PYA1bin(J1,1)
                        DUMMY2=-PSCALEFAC2vec(J2)!*PYA1bin(J2,1)
                   ELSE IF(k==3) THEN
                        DUMMY1=-PSCALEFAC2vec(J1)!*PYA1bin(J1,1)
                        DUMMY2=PSCALEFAC1vec(J2)!*PYA1bin(J2,1)
                   ELSE
                        DUMMY1=-PSCALEFAC2vec(J1)!*PYA1bin(J1,1)
                        DUMMY2=-PSCALEFAC2vec(J2)!*PYA1bin(J2,1)
                   END IF
                        ! first particle
                        xlj(k,1,:)=RI+DUMMY1*MATMUL(RMI,vecsbf)    ! vecsbf: (1,0,0) in the body frame of ellipsoid

                        ! second particle
                        xlj(k,2,:)=RJ+DUMMY2*MATMUL(RMJ,vecsbf)

                        ! separation between the LJ sites
                        rlj2(k)=(xlj(k,2,1)-xlj(k,1,1))**2+(xlj(k,2,2)-xlj(k,1,2))**2+(xlj(k,2,3)-xlj(k,1,3))**2
                        rlj(k)=sqrt(rlj2(k))
                        rljvec(k,1)=xlj(k,2,1)-xlj(k,1,1)
                        rljvec(k,2)=xlj(k,2,2)-xlj(k,1,2)
                        rljvec(k,3)=xlj(k,2,3)-xlj(k,1,3)

                        DUMMY=1.0D0/rlj(k)
                        rljunitvec(k,:)=rljvec(k,:)*DUMMY !/rlj(k)

                        drlj(k,1)=DUMMY*(xlj(k,2,1)-xlj(k,1,1))         !drlj/dx1
                        drlj(k,2)=DUMMY*(xlj(k,2,2)-xlj(k,1,2))         !drlj/dy1
                        drlj(k,3)=DUMMY*(xlj(k,2,3)-xlj(k,1,3))         !drlj/dz1
                        drlj(k,4)=-drlj(k,1)                               !drlj/dx2
                        drlj(k,5)=-drlj(k,2)                               !drlj/dy2
                        drlj(k,6)=-drlj(k,3)                               !drlj/dz2
                        drlj(k,7) =-DUMMY*DUMMY1*DOT_PRODUCT(rljvec(k,:),MATMUL(DPI1RM,vecsbf)) !drlj/dpx1
                        drlj(k,8) =-DUMMY*DUMMY1*DOT_PRODUCT(rljvec(k,:),MATMUL(DPI2RM,vecsbf)) !drlj/dpy1
                        drlj(k,9) =-DUMMY*DUMMY1*DOT_PRODUCT(rljvec(k,:),MATMUL(DPI3RM,vecsbf)) !drlj/dpz1
                        drlj(k,10) =  DUMMY*DUMMY2*DOT_PRODUCT(rljvec(k,:),MATMUL(DPJ1RM,vecsbf)) !drlj/dpx2
                        drlj(k,11) =  DUMMY*DUMMY2*DOT_PRODUCT(rljvec(k,:),MATMUL(DPJ2RM,vecsbf)) !drlj/dpy2
                        drlj(k,12) =  DUMMY*DUMMY2*DOT_PRODUCT(rljvec(k,:),MATMUL(DPJ3RM,vecsbf)) !drlj/dpz2

              ! interaction between the extra LJ sites:
                        LLJ(1,k)=sigma1(k)*DUMMY !/rlj(k)
                        LLJ(2,k)=LLJ(1,k)**2
                        LLJ(3,k)=LLJ(2,k)*LLJ(1,k)
                        LLJ(4,k)=LLJ(2,k)**2
                        LLJ(5,k)=LLJ(4,k)*LLJ(1,k)
                        LLJ(6,k)=LLJ(4,k)*LLJ(2,k)
                        LLJ(7,k)=LLJ(6,k)*LLJ(1,k)
                        LLJ(11,k)=LLJ(5,k)*LLJ(6,k)
                        LLJ(12,k)=LLJ(6,k)*LLJ(6,k)

!                            DUMMY=1.0D0/rlj(k)
!                            DUMMY=DUMMY**2
                            DO j=1,12
                                dLLJ1(k,j) =-sigma1(k)*DUMMY*DUMMY*drlj(k,j)
                            END DO

!                add corrections to
!                derivatives are zero at rc, and vanish smoothly with no discontinuities. 
!                VDUMMY=epsilon1*(LLJ12-0.0D0*LLJ6) !inner hard core
                 IF (PARAMONOVCUTOFF) THEN
                   IF (rlj(k)>=PCUTOFF) THEN
!                           VDUMMY = 0.0D0
                           term2(k)=1.0D0
                           term3(k)=0.0D0
                   ELSE

                        ! work out the spline terms for smooth cutoff of extra LJ sites
                        if(rlj2(k)<ron2) then
                           term2(k)=1.0D0
                           term3(k)=0.0D0
                        else if(rlj2(k)>ron2) then
                           term2(k)=(cut-rlj(k)**2)*(cut-rlj(k)**2)*(cut+2.0D0*rlj(k)**2-3.0D0*ron2)*range2inv3
                           term3(k)=rlj(k)*12.0D0*(cut-rlj(k)**2)*(ron2-rlj(k)**2)*range2inv3 ! d(term2)/dr
                        end if
                   END IF
!                   VDUMMY = 0.0D0
                 END IF ! IF (PARAMONOVCUTOFF)
!                 VDUMMY=4.0D0*epsilon1*(LLJ12-LLJ6)*term2 !extra LJ site
                 VDUMMY=VDUMMY+4.0D0*epsilon1(k,J1,J2)*term2(k)*(LLJ(12,k) - attr(k)*LLJ(6,k)) ! extra LJ sites (12-6)
!                 VDUMMY=VDUMMY+4.0D0*epsilon1(k,J1,J2)*term2(k)*LLJ(6,k) ! extra LJ sites (repulsive LJ6, testing)
              end do ! k=1,4
             END IF ! IF(LJSITE)

        IF(PARAMONOVCUTOFF) THEN
                !repulsive potential and periodic cutoff corrections
                VDUMMY = VDUMMY + 4.0D0 * PYEPSNOT * ( RHO112 + 6.0D0*CLJ12(1)*CLJ2(1)/RHO1SQ-7.0D0*CLJ12(1))
                    !                      LJ12(1)+(6.0D0*CLJ12(1))*CLJ2(1)/LJ2(1)-7.0D0*CLJ12(1)
                !attractive potential and periodic cutoff corrections
                VDUMMY = VDUMMY + 4.0D0 * PYEPSNOT * (- RHO26 - 3.0D0* CLJ6(2)*CLJ2(2)/RHO2SQ+4.0D0* CLJ6(2))
                    !                      (-LJ6(2)+(-3.0D0*CLJ6(2))*CLJ2(2)/LJ2(2) +4.0D0*CLJ6(2))
        ELSE
                VDUMMY = VDUMMY + 4.0D0 * PYEPSNOT * (RHO112 - 1.0D0 * RHO26)
        END IF

!       END IF ! IF(FROZEN)

               ENERGY = ENERGY + VDUMMY
               VT(J1) = VT(J1) + VDUMMY
               VT(J2) = VT(J2) + VDUMMY        ! pair potentials


            dVDUMMY(:) = 0.0D0

!     CALCULATE GRADIENT
             FIJ = 0.0D0
             TIJ = 0.0D0
             TJI = 0.0D0

 !       IF(.NOT.FROZEN(J1)) THEN
          IF(PARAMONOVCUTOFF) THEN
            dVDUMMY(:) = 0.0D0
            ! with respect to cartesians
               
            do j=1,3
              IF(LJSITE) THEN
               do k=1,maxinteractions
!               dVDUMMY(j) = 24.0D0*epsilon1*(2.0D0*LLJ11*dLLJ1(j)-LLJ5*dLLJ1(j))*1.0D0*term2 + &
!                          & 4.0D0*epsilon1*(LLJ12-LLJ6)*term3*drlj(j)! extra LJ site derivatives
                dVDUMMY(j) = dVDUMMY(j) + 4.0D0*epsilon1(k,J1,J2)*(12.0D0*LLJ(11,k)*dLLJ1(k,j))*term2(k) + &
                          & 4.0D0*epsilon1(k,J1,J2)*LLJ(12,k)*term3(k)*drlj(k,j)! extra LJ site derivatives, now only repulsive
                dVDUMMY(j) = dVDUMMY(j) - attr(k)*(4.0D0*epsilon1(k,J1,J2)*(6.0D0*LLJ(5,k)*dLLJ1(k,j))*term2(k) + &
                          & 4.0D0*epsilon1(k,J1,J2)*LLJ(6,k)*term3(k)*drlj(k,j)) ! attractive secondary apex site
               end do
              ELSE
               dVDUMMY(j) = 0.0D0
              END IF
!                dVDUMMY(j) = dLLJ1(j)
               dVDUMMY(j) = dVDUMMY(j) + 24.0D0 * PYEPSNOT * ( &
               & 2.0D0*RHO1SQ*DG1DR(j)*(RHO16*RHO1SQ*RHO1SQ*RHO1-CLJ14(1)/(RHO1SQ*RHO1))&
               &+14.0D0*dCLJ1(1,j)*(CLJ13(1)/RHO1SQ-CLJ11(1))) 
               !repulsive derivative

               dVDUMMY(j) = dVDUMMY(j) + 24.0D0 * PYEPSNOT * (-1.0D0*RHO2SQ*DG2DR(j)*(RHO2SQ*RHO2SQ*RHO2-CLJ8(2)/(RHO2SQ*RHO2))&
                   &-4.0D0*dCLJ1(2,j)*(CLJ7(2)/RHO2SQ-CLJ5(2))) !attractive derivative
!               dVDUMMY(j) = dCLJ1(1,j)
               FIJ(j) = dVDUMMY(j)
            end do

            do k=1,2
                DUMMY=LJ1(k)/PYSIGNOT
                DUMMY=DUMMY**2
                DUMMY1=SRTFI(k)
                DUMMY2=DUMMY1**3
               do j=7,12
                  dLJ1(k,j) =-1.0D0*PYSIGNOT*DUMMY*(0.5D0*ABSRIJ*DUMMY2*DFP(k,j-6))
               end do
            end do
            ! with respect to orientation vectors
            do j=7,12
             IF(LJSITE) THEN
              do k=1,maxinteractions
!               dVDUMMY(j) = 24.0D0*epsilon1*(2.0D0*LLJ1**11*dLLJ1(j)-LLJ1**5*dLLJ1(j))*term2 + &
!                               & 4.0D0*epsilon1*(LLJ12-LLJ6)*term3*drlj(j)! extra LJ site derivatives
               dVDUMMY(j) = dVDUMMY(j) + 4.0D0*epsilon1(k,J1,J2)*(12.0D0*LLJ(11,k)*dLLJ1(k,j))*term2(k) + &
                               & 4.0D0*epsilon1(k,J1,J2)*LLJ(12,k)*term3(k)*drlj(k,j)! extra LJ site derivatives, now only repulsive
               dVDUMMY(j) = dVDUMMY(j) - attr(k)*(4.0D0*epsilon1(k,J1,J2)*(6.0D0*LLJ(5,k)*dLLJ1(k,j))*term2(k) + &
                          & 4.0D0*epsilon1(k,J1,J2)*LLJ(6,k)*term3(k)*drlj(k,j)) ! attractive secondary apex site
              end do
             ELSE
               dVDUMMY(j) = 0.0
             END IF
!                dVDUMMY(j) = dLLJ1(j)
               dVDUMMY(j) = dVDUMMY(j) + 24.0D0 * PYEPSNOT * ( 2.0D0*dLJ1(1,j)*(RHO16*RHO1SQ*RHO1SQ*RHO1-CLJ14(1)/(RHO1SQ*RHO1))&
                   &+14.0D0*dCLJ1(1,j)*(CLJ13(1)/RHO1SQ-CLJ11(1))) !repulsive derivative

               dVDUMMY(j) = dVDUMMY(j) + 24.0D0 * PYEPSNOT * (-1.0D0*dLJ1(2,j)*(RHO2SQ*RHO2SQ*RHO2- CLJ8(2)/(RHO2SQ*RHO2))&
                   &- 4.0D0*dCLJ1(2,j)*( CLJ7(2)/RHO2SQ- CLJ5(2))) !attractive derivative
!               write(*,*) 'dLJ1(2,j)', DFP(2,j-6),j-6
!               dVDUMMY(j) = dCLJ1(1,j)
            end do
            do j=1,3
               TIJ(j) = dVDUMMY(6+j)
               TJI(j) = dVDUMMY(9+j)
            end do
!        ELSE
                 
!        DO j=1,12
!          dvdummy(j)=drlj(j)
!        END DO

          ELSE  !no cutoff
             FIJ = 0.0D0
             TIJ = 0.0D0
             TJI = 0.0D0
             dVDUMMY(:)=0.0D0

           IF(LJSITE) THEN
            do k=1,maxinteractions
             do j=1,3
               dVDUMMY(j) = dVDUMMY(j) + 4.0D0*epsilon1(k,J1,J2)*(12.0D0*LLJ(11,k)*dLLJ1(k,j))*term2(k) + &
                          & 4.0D0*epsilon1(k,J1,J2)*LLJ(12,k)*term3(k)*drlj(k,j)! extra LJ site derivatives, currently only repulsive (LJ12)
!               dVDUMMY(j) = dVDUMMY(j) + 4.0D0*epsilon1(k,J1,J2)*(6.0D0*LLJ(5,k)*dLLJ1(k,j))*term2(k)! + &
!!                          & 4.0D0*epsilon1(k,J1,J2)*LLJ(12,k)*term3(k)*drlj(k,j)! extra LJ site derivatives, currently only repulsive (LJ6)
               dVDUMMY(j) = dVDUMMY(j) - attr(k)*(4.0D0*epsilon1(k,J1,J2)*(6.0D0*LLJ(5,k)*dLLJ1(k,j))*term2(k) + &
                          & 4.0D0*epsilon1(k,J1,J2)*LLJ(6,k)*term3(k)*drlj(k,j)) ! attractive secondary apex site
               FIJ(j) = dVDUMMY(j)
             end do
             do j=7,12
               dVDUMMY(j) = dVDUMMY(j) + 4.0D0*epsilon1(k,J1,J2)*(12.0D0*LLJ(11,k)*dLLJ1(k,j))*term2(k) + &
                          & 4.0D0*epsilon1(k,J1,J2)*LLJ(12,k)*term3(k)*drlj(k,j)! extra LJ site derivatives, currently only repulsive
               dVDUMMY(j) = dVDUMMY(j) - attr(k)*(4.0D0*epsilon1(k,J1,J2)*(6.0D0*LLJ(5,k)*dLLJ1(k,j))*term2(k) + &
                          & 4.0D0*epsilon1(k,J1,J2)*LLJ(6,k)*term3(k)*drlj(k,j)) ! attractive secondary apex site
             end do
            end do
           END IF !if LJSITE 
!        END IF !IF(.NOT.FROZEN)
!4.0D0*epsilon0*(12.0D0*LJ11(1)*dLJ1(1,j)-6.0D0*LJ5(2)*dLJ1(2,j))
               FIJ    = FIJ + 24.0D0*PYEPSNOT*(2.D0*RHO112*RHO1*DG1DR - 1.0D0 * RHO26*RHO2*DG2DR)
               TIJ(1) = DVDF1*DF1PI1 + 1.0D0* DVDF2*DF2PI1
               TIJ(2) = DVDF1*DF1PI2 + 1.0D0*  DVDF2*DF2PI2
               TIJ(3) = DVDF1*DF1PI3 + 1.0D0*  DVDF2*DF2PI3
               TJI(1) = DVDF1*DF1PJ1 + 1.0D0*  DVDF2*DF2PJ1
               TJI(2) = DVDF1*DF1PJ2 + 1.0D0*  DVDF2*DF2PJ2
               TJI(3) = DVDF1*DF1PJ3 + 1.0D0*  DVDF2*DF2PJ3
               TIJ = dVDUMMY(7:9) + 24.0D0 * PYEPSNOT * TIJ
               TJI = dVDUMMY(10:12) + 24.0D0 * PYEPSNOT * TJI
          END IF ! cutoff/no cutoff block

!                G(J3-2)=G(J3-2)+dVDUMMY(1)
!                G(j4-2)=G(j4-2)+dVDUMMY(4)
!                G(J3-1)=G(J3-1)+dVDUMMY(2)
!                G(j4-1)=G(j4-1)+dVDUMMY(5)
!                G(J3  )=G(J3  )+dVDUMMY(3)
!                G(j4  )=G(j4  )+dVDUMMY(6)



!              G(J3-2:J3) = G(J3-2:J3) + dvdummy(1:3)
!              G(J4-2:J4) = G(J4-2:J4) + dvdummy(4:6)
!              G(J5-2:J5) = G(J5-2:J5) + dvdummy(7:9)
!              G(J6-2:J6) = G(J6-2:J6) + dvdummy(10:12)

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



! PY potential, sf344's implementation
       SUBROUTINE OLDPARAMONOV(X,V,EGB,GTEST,STEST)
       use commons, only: BOXLX,BOXLY,BOXLZ,PARAMONOVPBCX,PARAMONOVPBCY,PARAMONOVPBCZ,PCUTOFF,PARAMONOVCUTOFF,&
                &       natoms,psigma0,pepsilon0,parama1,paramb1,paramc1,parama2,paramb2,paramc2,VT,LJSITE,PEPSILON1
       !use keyword
       IMPLICIT NONE
!      Anisotropic potential by Paramonov et. al., JCP 123, 194111 (2005)
!      finally working for triaxial systems as well

       DOUBLE PRECISION :: x1,y1,z1,x2,y2,z2,a1(2),b1(2),c1(2),F1(2),dF1(2,12),r,rlj,rlj2,xlj(2,3),drlj(12), &
                        &  term2, term3,ron2,ron,cut,range2inv3
       DOUBLE PRECISION :: px1(2),py1(2),pz1(2),alpha(2),beta(2),gamma(2), &
                        &                   sa(2),sb(2),sg(2),ca(2),cb(2),cg(2),&
                        &                   px(2),py(2),pz(2),qx(2),qy(2),qz(2),rx(2),ry(2),rz(2)
       DOUBLE PRECISION :: da(2,3),db(2,3),dg(2,3),dpx(2,3),dpy(2,3),dpz(2,3),dqx(2,3),dqy(2,3),&
                        &                   dqz(2,3),drx(2,3),dry(2,3),drz(2,3)
       DOUBLE PRECISION :: LJ1(2),LJ2(2),LJ3(2),LJ5(2),LJ6(2), &
                        &                  LJ12(2), LJ4(2), LJ7(2), LJ11(2),dLJ1(2,12),dr(6)
       DOUBLE PRECISION :: lambda(2), A11(3,3),B11(3,3),C_1(6),C1inv(3,3),W(3),DUMMY,stepsize,dF_lambdaold,lambdaold
       DOUBLE PRECISION :: DUMMY1, DUMMY2, S_DUMMY1, C_DUMMY1, S_DUMMY2, C_DUMMY2, VDUMMY, dVDUMMY(12),epsilon0
       DOUBLE PRECISION :: DUMMYX, DUMMYY,DUMMYZ, sigma0(2), dF_lambda(1), EGB,X(3*NATOMS), &
                        &  V(3*NATOMS),DPRAND,convg,u1(3),u2(3),R1(3),R1T(1,3)
       INTEGER   :: i,j,k,l,iter=0,REALNATOMS,J1,J2,J3,K2,J4
       LOGICAL   :: GTEST,STEST
       DOUBLE PRECISION :: G_invR(3),G_invRT(1,3),S_lambda(1), RANDOM,PI,Fcheck(3),Dcheck
       DOUBLE PRECISION :: dGu1(3,3,3),dGu2(3,3,3),dXu1(3,3),dXu2(3,3),dFu1(3,1),dFu2(3,1)
       DOUBLE PRECISION :: dAA(3,2,3,3),dBB(3,2,3,3)

       DOUBLE PRECISION :: LLJ1,LLJ2,LLJ3,LLJ4,LLJ5,LLJ6,LLJ7,LLJ11,LLJ12,dLLJ1(12),sigma1,epsilon1
       DOUBLE PRECISION :: CLJ1(2),CLJ2(2),CLJ3(2),CLJ4(2),CLJ5(2),CLJ6(2)
       DOUBLE PRECISION :: CLJ7(2),CLJ8(2),CLJ11(2),CLJ12(2),CLJ13(2),CLJ14(2),dCLJ1(2,12)

       epsilon1=PEPSILON1(1)   ! we are going to use sigma1 and epsilon1 for the extra LJ sites
       sigma1=1.0D0!PSIGMA1 is nonexistent from now on
       REALNATOMS=NATOMS/2
       PI=4.0D0*ATAN(1.0)
! cutoff stuff for the extra LJ sites
       cut = PCUTOFF*PCUTOFF ! cutoff squared
       ron = PCUTOFF*0.9D0
       ron2 = ron*ron
       range2inv3=1.0D0/(cut-ron2)**3
       term2=1.0D0
       term3=0.0D0
!      VT(1:REALNATOMS)=0.0D0
       V(1:6*REALNATOMS)=0.0D0
!      WRITE(*,*) "In energy cycle"
!      WRITE(*,*) X(:)
!      write(*,*) BOXLX, BOXLY, BOXLZ, PCUTOFF
!      write (*,*) PARAMONOVPBCX,PARAMONOVPBCY,PARAMONOVPBCZ, PARAMONOVCUTOFF

       EGB=0.0D0

       sigma0(1)=PSIGMA0(1)
       sigma0(2)=PSIGMA0(2)
       epsilon0=PEPSILON0

       a1(1)=PARAMa1**2
       b1(1)=PARAMb1**2
       c1(1)=PARAMc1**2

       a1(2)=PARAMa2**2
       b1(2)=PARAMb2**2
       c1(2)=PARAMc2**2

!      WRITE(*,*) a1(:),b1(:),c1(:)

       DO J1=1,REALNATOMS
        J2=3*J1

        IF (PARAMONOVPBCX) THEN
                !ensure x component of particle 1 vector is within BoxLx/2 of zero. If it isn't then subtract integer number of boxlx's such that it is.
                X(J2-2)=X(J2-2)-BOXLX*NINT(X(J2-2)/BOXLX)
        ENDIF

        IF (PARAMONOVPBCY) THEN
                !ensure y component of particle 1 vector is within BoxLy/2 of zero. If it isn't then subtract integer number of boxly's such that it is.
                X(J2-1)=X(J2-1)-BOXLY*NINT(X(J2-1)/BOXLY)
        END IF

        IF (PARAMONOVPBCZ) THEN
                !ensure z component of particle 1 vector is within BoxLz/2 of zero. If it isn't then subtract integer number of boxlz's such that it is.
                X(J2)=X(J2)-BOXLZ*NINT(X(J2)/BOXLZ)
        ENDIF


         
       END DO

       DO J1=1,REALNATOMS
        J2=3*J1
        
        x1 = X(J2-2)
        y1 = X(J2-1)
        z1 = X(J2  )

!        px, py, pz: angle-axis coordinate system (BUT: right hand rotation!)
        px1(1) = X(3*REALNATOMS+J2-2)
        py1(1) = X(3*REALNATOMS+J2-1)
        pz1(1) = X(3*REALNATOMS+J2  )
!        IF(px1(1)==0.0D0) px1(1)=px1(1)+4.0D0*pi**2
!        IF(py1(1)==0.0D0) py1(1)=py1(1)+4.0D0*pi**2
!        IF(pz1(1)==0.0D0) pz1(1)=pz1(1)+4.0D0*pi**2

        DO J3=J1+1,REALNATOMS
                VDUMMY=0.0D0
                K2=3*J3
!                IF (PARAMONOVPBCX) THEN
!                   !ensure x component of particle 2 vector is within BoxLx/2 of zero. If it isn't then subtract integer number of boxlx's such that it is.
!                        X(K2-2)=X(K2-2)-BOXLX*NINT(X(K2-2)/BOXLX)
!                ENDIF
!
!                IF (PARAMONOVPBCY) THEN
!                   !ensure y component of particle 2 vector is within BoxLy/2 of zero. If it isn't then subtract integer number of boxly's such that it is.
!                        X(K2-1)=X(K2-1)-BOXLY*NINT(X(K2-1)/BOXLY)
!                ENDIF
!
!                IF (PARAMONOVPBCZ) THEN
!                   !ensure z component of particle 2 vector is within BoxLz/2 of zero. If it isn't then subtract integer number of boxlz's such that it is.
!                        X(K2)=X(K2)-BOXLZ*NINT(X(K2)/BOXLZ)
!                ENDIF
                x2 = X(K2-2)
                y2 = X(K2-1)
                z2 = X(K2  )
                px1(2) = X(3*REALNATOMS+K2-2)
                py1(2) = X(3*REALNATOMS+K2-1)
                pz1(2) = X(3*REALNATOMS+K2  )
!                IF(px1(2)==0.0D0) px1(2)=px1(2)+4.0D0*pi**2
!                IF(py1(2)==0.0D0) py1(2)=py1(2)+4.0D0*pi**2
!                IF(pz1(2)==0.0D0) pz1(2)=pz1(2)+4.0D0*pi**2

!            derivative check
!             DO l=1,3
 
                R1(1)=x2-x1
                R1(2)=y2-y1
                R1(3)=z2-z1

                !compute vector to image of second particle that is within boxlx/2 of first particle
                IF (PARAMONOVPBCX) THEN
                        !write(*,*) "x1,x2,r1(1)"
                        !write(*,*) x1,x2,R1(1)
                        !find the x-component of the vector to the nearest image of the second particle
                        R1(1)=R1(1)-BOXLX*NINT(R1(1)/BOXLX)
                        
                        !recompute x2 with the new component. Don't worry if it's outside the box.
                        x2=R1(1)+x1
                        !write(*,*) "x1,x2,r1(1)"
                        !write(*,*) x1,x2,R1(1)
                ENDIF

                !Set periodic boundary conditions for Y dimension.
                IF (PARAMONOVPBCY) THEN
                        !write(*,*) "y1,y2,r1(2)"
                        !write(*,*) y1,y2,R1(2)

                        !find the y-component of the vector to the nearest image of the second particle
                        R1(2)=R1(2)-BOXLY*NINT(R1(2)/BOXLY)

                        !recompute y2 with the new component. Don't worry if it's outside the box.
                        y2=R1(2)+y1
                        !write(*,*) "y1,y2,r1(2)"
                        !write(*,*) y1,y2,R1(2)
                ENDIF

                !Set periodic boundary conditions for the z dimension
                IF (PARAMONOVPBCZ) THEN
                        !write(*,*) "z1,z2,r1(3)"
                        !write(*,*) z1,z2,R1(3)
                        !Find the z-component of the vector to the nearest image of the second particle
                        R1(3)=R1(3)-BOXLZ*NINT(R1(3)/BOXLZ)
                        !Recompute z2 with the new component.
                        z2=R1(3)+z1
                        !write(*,*) "z1,z2,r1(3)"
                        !write(*,*) z1,z2,R1(3)
                END IF

                !check for cold fusion
                IF(SQRT((x2-x1)**2+(y2-y1)**2+(z2-z1)**2).LE.0.01D0) THEN
                        WRITE(*,*) 'sf344> cold fusion detected'
                        VDUMMY=-1.0D20
                        EGB=-1.0D20
                        V(:)=0.0D0
                        GOTO 999
                END IF
                                
                r=SQRT((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)

!if we are cutting off the potential and the inter particle distance is >=PCUTOFF don't bother doing any maths. the answer is zero.
                IF ((.NOT.PARAMONOVCUTOFF).OR.(PARAMONOVCUTOFF.AND.(r<PCUTOFF))) THEN
                        R1T(1,1)=R1(1)
                        R1T(1,2)=R1(2)
                        R1T(1,3)=R1(3)

                        dr(1)=-1.0D0/r*(x2-x1)
                        dr(2)=-1.0D0/r*(y2-y1)
                        dr(3)=-1.0D0/r*(z2-z1)
                        dr(4)=-dr(1)
                        dr(5)=-dr(2)
                        dr(6)=-dr(3)
                          
                        DO k=1,2

    
                                stepsize=1.0D0
                                lambda(k)=0.5D0
                                IF (k==2) THEN
                                        lambda(2)=lambda(1)
                                END IF
                                ITER=0
                                convg=1.0D-8
                                dF_lambdaold=1.0D6
 
                                DO i=1,2
                                        gamma(i) = sqrt(px1(i)**2+py1(i)**2+pz1(i)**2)
!                                        write(*,*) 'gamma(i) before',i,gamma(i)
!                                        gamma(i) = gamma(i)-2.0D0*PI*NINT(gamma(i)/(2.0D0*PI))+2.0D0*PI
!                                        write(*,*) 'gamma(i) after',i,gamma(i)

                                        dg(i,1) = px1(i)/gamma(i)
                                        dg(i,2) = py1(i)/gamma(i)
                                        dg(i,3) = pz1(i)/gamma(i)
                                        px(i) = px1(i)/gamma(i)
                                        py(i) = py1(i)/gamma(i)
                                        pz(i) = pz1(i)/gamma(i)
alpha(i)=ATAN2(py(i),px(i))
!                                                                da(i,1) = -1.0D0*py1(i)/(px1(i)**2+py1(i)**2)
!                                                                da(i,2) = px1(i)/(px1(i)**2+py1(i)**2)

                                                                da(i,1) = -1.0D0*py1(i)/(px1(i)**2+py1(i)**2)
                                                                da(i,2) = px1(i)/(gamma(i)**2-pz1(i)**2)
                                                                da(i,3) = 0.0d0
 
!IF(0) THEN
                                        IF(px(i).eq.0.0d0) THEN
                                                IF(py(i)>=0.0D0) THEN
                                                        alpha(i) = PI/2           ! Euler angle alpha 
                                                ELSE
                                                        alpha(i) = -1.0D0*PI/2
                                                END IF
                                        ELSE
                                                IF(py(i)>=0.0D0) THEN
                                                        IF(px(i)>0.0D0) THEN
                                                                alpha(i) = 1.0D0*atan(py(i)/px(i))   ! first quadrant
                                                                da(i,1) = -1.0D0*py1(i)/(px1(i)**2+py1(i)**2)
                                                                da(i,2) = px1(i)/(px1(i)**2+py1(i)**2)
                                                                da(i,3) = 0.0d0
                                                        ELSE    ! px<0
                                                                alpha(i) = 1.0D0*atan(py(i)/px(i))+PI ! should be in the second quadrant
                                                                da(i,1) = -1.0D0*py1(i)/(px1(i)**2+py1(i)**2)
                                                                da(i,2) = px1(i)/(px1(i)**2+py1(i)**2)
                                                                da(i,3) = 0.0d0
                                                        END IF
                                                ELSE IF(py(i)<0.0D0) THEN
                                                        IF(px(i)>0.0D0) THEN
                                                                alpha(i) = 1.0D0*atan(py(i)/px(i))    ! fourth quadrant
                                                                da(i,1) = -1.0D0*py1(i)/(px1(i)**2+py1(i)**2)
                                                                da(i,2) = px1(i)/(px1(i)**2+py1(i)**2)
                                                                da(i,3) = 0.0d0
                                                        ELSE    ! px<0
                                                                alpha(i) = 1.0D0*atan(py(i)/px(i))-PI            ! third quadrant
                                                                da(i,1) = -1.0D0*py1(i)/(px1(i)**2+py1(i)**2)
                                                                da(i,2) = px1(i)/(px1(i)**2+py1(i)**2)
                                                                da(i,3) = 0.0d0
                                                        END IF
                                                END IF 
                                        END IF
!END IF 
                                        beta(i) = 1.0D0*acos(pz(i))
                                                DUMMY1 = pz1(i)/(sqrt(px1(i)**2+py1(i)**2)*gamma(i)**2)
                                        db(i,1) = px1(i)*DUMMY1
                                        db(i,2) = py1(i)*DUMMY1
                                        db(i,3) = -1.0D0*sqrt(px1(i)**2+py1(i)**2)/gamma(i)**2
        
                                        sa(i) = sin(alpha(i))
                                        sb(i) = sin(beta(i))
                                        sg(i) = sin(gamma(i))
                                        ca(i) = cos(alpha(i))
                                        cb(i) = cos(beta(i))
                                        cg(i) = cos(gamma(i))

                                        qx(i) = -sg(i)*cb(i)*ca(i)-cg(i)*sa(i)
                                        qy(i) = cg(i)*ca(i)-sg(i)*cb(i)*sa(i)
                                        qz(i) = sg(i)*sb(i)
                                        rx(i) = cg(i)*cb(i)*ca(i)-sg(i)*sa(i)
                                        ry(i) = cg(i)*cb(i)*sa(i)+sg(i)*ca(i)
                                        rz(i) = -cg(i)*sb(i)


                                        DO j=1,3
                                                dqx(i,j) = -cg(i)*(cb(i)*ca(i))*dg(i,j)-sg(i)*(-sb(i)*ca(i)*db(i,j)-&
                                                          & cb(i)*sa(i)*da(i,j))-(-sg(i)*dg(i,j)*sa(i)+cg(i)*ca(i)*da(i,j))
                                                dqy(i,j) = -sg(i)*dg(i,j)*ca(i)-cg(i)*sa(i)*da(i,j)-(cg(i)*cb(i)*sa(i)*dg(i,j)+&
                                                          & sg(i)*(-sb(i)*db(i,j)*sa(i)+cb(i)*ca(i)*da(i,j)))
                                                dqz(i,j) = cg(i)*sb(i)*dg(i,j)+sg(i)*cb(i)*db(i,j)
                  drx(i,j) = -sg(i)*dg(i,j)*cb(i)*ca(i)+cg(i)*(-sb(i)*ca(i)*db(i,j)-cb(i)*sa(i)*da(i,j))-&
                                                          & (cg(i)*sa(i)*dg(i,j)+sg(i)*ca(i)*da(i,j))
                  dry(i,j) = -sg(i)*dg(i,j)*cb(i)*sa(i)+cg(i)*(-sb(i)*sa(i)*db(i,j)+cb(i)*ca(i)*da(i,j))+&
                                                          &  cg(i)*ca(i)*dg(i,j)-sg(i)*sa(i)*da(i,j)
                                                drz(i,j) = sg(i)*sb(i)*dg(i,j)-cg(i)*cb(i)*db(i,j)
                                        END DO ! j=1,3
                                END DO  ! i=1,2

! extra LJ site at the end of the repulsive semiaxis C (unit vector pointing to it: (rx,ry,rz))
                                        IF(LJSITE.AND.k==1) THEN
                                             ! first particle
                                                xlj(1,1)=x1+PARAMc1*rx(1)
                                                xlj(1,2)=y1+PARAMc1*ry(1)
                                                xlj(1,3)=z1+PARAMc1*rz(1)
                                             ! second particle
                                                xlj(2,1)=x2+PARAMc1*rx(2)
                                                xlj(2,2)=y2+PARAMc1*ry(2)
                                                xlj(2,3)=z2+PARAMc1*rz(2)
                                             DO i=1,3
                                                rlj2=(xlj(2,1)-xlj(1,1))**2+(xlj(2,2)-xlj(1,2))**2+(xlj(2,3)-xlj(1,3))**2
                                                rlj=sqrt(rlj2)
                                             END DO
                                                drlj(1)=-1.0D0/rlj*(xlj(2,1)-xlj(1,1))         !drlj/dx1
                                                drlj(2)=-1.0D0/rlj*(xlj(2,2)-xlj(1,2))         !drlj/dy1
                                                drlj(3)=-1.0D0/rlj*(xlj(2,3)-xlj(1,3))         !drlj/dz1
                                                drlj(4)=-drlj(1)                               !drlj/dx2
                                                drlj(5)=-drlj(2)                               !drlj/dy2
                                                drlj(6)=-drlj(3)                               !drlj/dz2

                                                drlj(7) = 1.0D0/rlj*( (xlj(2,1)-xlj(1,1))*PARAMc1*(-drx(1,1)) + &
                                                         &            (xlj(2,2)-xlj(1,2))*PARAMc1*(-dry(1,1)) + &
                                                         &            (xlj(2,3)-xlj(1,3))*PARAMc1*(-drz(1,1))   ) !drlj/dpx1
                                                drlj(8) = 1.0D0/rlj*( (xlj(2,1)-xlj(1,1))*PARAMc1*(-drx(1,2)) + &
                                                         &            (xlj(2,2)-xlj(1,2))*PARAMc1*(-dry(1,2)) + &
                                                         &            (xlj(2,3)-xlj(1,3))*PARAMc1*(-drz(1,2))   ) !drlj/dpy1
                                                drlj(9) = 1.0D0/rlj*( (xlj(2,1)-xlj(1,1))*PARAMc1*(-drx(1,3)) + &
                                                         &            (xlj(2,2)-xlj(1,2))*PARAMc1*(-dry(1,3)) + &
                                                         &            (xlj(2,3)-xlj(1,3))*PARAMc1*(-drz(1,3))   ) !drlj/dpz1
                                                drlj(10)= 1.0D0/rlj*( (xlj(2,1)-xlj(1,1))*PARAMc1*( drx(2,1)) + &
                                                         &            (xlj(2,2)-xlj(1,2))*PARAMc1*( dry(2,1)) + &
                                                         &            (xlj(2,3)-xlj(1,3))*PARAMc1*( drz(2,1))   ) !drlj/dpx2
                                                drlj(11)= 1.0D0/rlj*( (xlj(2,1)-xlj(1,1))*PARAMc1*( drx(2,2)) + &
                                                         &            (xlj(2,2)-xlj(1,2))*PARAMc1*( dry(2,2)) + &
                                                         &            (xlj(2,3)-xlj(1,3))*PARAMc1*( drz(2,2))   ) !drlj/dpy2
                                                drlj(12)= 1.0D0/rlj*( (xlj(2,1)-xlj(1,1))*PARAMc1*( drx(2,3)) + &
                                                         &            (xlj(2,2)-xlj(1,2))*PARAMc1*( dry(2,3)) + &
                                                         &            (xlj(2,3)-xlj(1,3))*PARAMc1*( drz(2,3))   ) !drlj/dpz2

                                        END IF

                        A11(1,1)=px(1)**2*a1(k)+qx(1)**2*b1(k)+rx(1)**2*c1(k)
                        A11(1,2)=py(1)*px(1)*a1(k)+qy(1)*qx(1)*b1(k)+ry(1)*rx(1)*c1(k)
                        A11(1,3)=px(1)*pz(1)*a1(k)+qx(1)*qz(1)*b1(k)+rx(1)*rz(1)*c1(k)
                        A11(2,1)=A11(1,2)
                        A11(2,2)=py(1)**2*a1(k)+qy(1)**2*b1(k)+ry(1)**2*c1(k)
                        A11(2,3)=py(1)*pz(1)*a1(k)+qy(1)*qz(1)*b1(k)+ry(1)*rz(1)*c1(k)
                        A11(3,1)=A11(1,3)
                        A11(3,2)=A11(2,3)
                        A11(3,3)=pz(1)**2*a1(k)+qz(1)**2*b1(k)+rz(1)**2*c1(k)

                        B11(1,1)=px(2)**2*a1(k)+qx(2)**2*b1(k)+rx(2)**2*c1(k)
                        B11(1,2)=py(2)*px(2)*a1(k)+qy(2)*qx(2)*b1(k)+ry(2)*rx(2)*c1(k)
                        B11(1,3)=px(2)*pz(2)*a1(k)+qx(2)*qz(2)*b1(k)+rx(2)*rz(2)*c1(k)
                        B11(2,1)=B11(1,2)
                        B11(2,2)=py(2)**2*a1(k)+qy(2)**2*b1(k)+ry(2)**2*c1(k)
                        B11(2,3)=py(2)*pz(2)*a1(k)+qy(2)*qz(2)*b1(k)+ry(2)*rz(2)*c1(k)
                        B11(3,1)=B11(1,3)
                        B11(3,2)=B11(2,3)
                        B11(3,3)=pz(2)**2*a1(k)+qz(2)**2*b1(k)+rz(2)**2*c1(k)

                        DO
                                RANDOM=DPRAND()
                                 dF_lambda(:)=0.0D0
                                       W(:)=0.0D0

                                C_1(1)=(1.0D0-lambda(k))*A11(1,1)+lambda(k)*B11(1,1)
                                C_1(2)=(1.0D0-lambda(k))*A11(1,2)+lambda(k)*B11(1,2)
                                C_1(3)=(1.0D0-lambda(k))*A11(1,3)+lambda(k)*B11(1,3)
                                C_1(4)=(1.0D0-lambda(k))*A11(2,2)+lambda(k)*B11(2,2)
                                        C_1(5)=(1.0D0-lambda(k))*A11(2,3)+lambda(k)*B11(2,3)
                                        C_1(6)=(1.0D0-lambda(k))*A11(3,3)+lambda(k)*B11(3,3)

                                CALL SVERT(C_1,3,W)  !inverting the symmetric matrix C1 (or G, eq.2.6 in ECF paper)
     
                                C1inv(1,1)=C_1(1)
                                C1inv(1,2)=C_1(2)
                                C1inv(1,3)=C_1(3)
                                C1inv(2,1)=C1inv(1,2)
                                C1inv(2,2)=C_1(4)
                                C1inv(2,3)=C_1(5)
                                C1inv(3,1)=C1inv(1,3)
                                C1inv(3,2)=C1inv(2,3)
                                C1inv(3,3)=C_1(6)

                                G_invR=MATMUL(C1inv,R1)
     
                                G_invRT(1,1)=G_invR(1)
                                G_invRT(1,2)=G_invR(2)
                                G_invRT(1,3)=G_invR(3)
     
                                dF_lambda=MATMUL(G_invRT,MATMUL((1.0D0-lambda(k))**2*A11-lambda(k)**2*B11,G_invR)) !derivative of the contact function with respect to lambda
                                DUMMY=dF_lambda(1)
                                DUMMY1=lambda(k)

                                IF(ABS(dF_lambda(1))<convg) THEN
                                        EXIT
                                ELSE
                                        ITER=ITER+1
!                                        if the new derivative is larger than the old one, reset lambda to its previous value and decrease stepsize
                                        IF(ABS(dF_lambda(1))>ABS(dF_lambdaold)) THEN
!                                                WRITE(*,*) 'dF_lambda>df_lambdaold', dF_lambda, df_lambdaold
                                                lambda(k)=lambdaold
                                                dF_lambda(1)=dF_lambdaold
!                                                IF(stepsize>=0.000000001D0)
                                                stepsize=stepsize/3.0D0
                                        END IF

                                        IF(lambda(k)<0.OR.lambda(k)>1) THEN
!                                                WRITE(*,*) "lambda value not in [0,1], resetting",J1,J3
                                                stepsize=1.0D0
                                                lambda(k)=0.5
                                                dF_lambda(1)=0.1D0
                                                DUMMY=dF_lambda(1)
                                                DUMMY1=lambda(k)
                                        END IF

                                        IF(ABS(dF_lambda(1))<=0.1D0) THEN
                                                lambda(k)=lambda(k)+dF_lambda(1)*stepsize*RANDOM
                                        ELSE IF(ABS(dF_lambda(1))>1.0D0) THEN
                                                        lambda(k)=lambda(k)+(dF_lambda(1)/ABS(dF_lambda(1)))*RANDOM*stepsize*0.2
                                        ELSE
                                                        lambda(k)=lambda(k)+dF_lambda(1)*stepsize*0.1*RANDOM
                                        END IF

                                        lambdaold=DUMMY1
                                        dF_lambdaold=DUMMY
 
                                END IF
!                                IF(ITER.GT.1000) THEN
!                                        WRITE(*,*) lambda(k),dF_lambda(1),ITER,px,py,py/px,alpha
!                                        WRITE(*,*) x1,y1,z1
!                                        WRITE(*,*) x2,y2,z2
!                                        WRITE(*,*) px1,py1,pz1
!                                        WRITE(*,*) px2,py2,pz2
!                                        WRITE(*,*) SQRT(px1**2+py1**2+pz1**2)
!                                        WRITE(*,*) SQRT(px2**2+py2**2+pz2**2)
!                                END IF
                        END DO

! calculate S(lambda) at the cutoff, assuming lambda=lambda_c

                        S_lambda=lambda(k)*(1.0D0-lambda(k))*MATMUL(R1T,G_invR)  ! this is the ellipsoid contact function
                        F1(k)=S_lambda(1)
!                        DERIVATIVES of the ELLIPSOID CONTACT FUNCTION

!                        with respect to cartesian coordinates:
                        DUMMY=2.0D0*lambda(k)*(1.0D0-lambda(k))
                            DO j=1,3
                                dF1(k,j)=DUMMY*G_invR(j)*(-1.0D0)  !dF/dx1,dy1,dz1
                                dF1(k,j+3)=-dF1(k,j)           !dF/dx2,dy2,dz2
                            END DO

!                        with respect to orientation:
                        DO j=1,2
                                dpx(j,1) = (py1(j)**2+pz1(j)**2)/gamma(j)**3                      ! dpx/dpx1
                                dpx(j,2) = -1.0D0*(px1(j)*py1(j))/gamma(j)**3                     ! dpx/dpy1
                                dpx(j,3) = -1.0D0*(px1(j)*pz1(j))/gamma(j)**3                     ! dpx/dpz1
                                dpy(j,1) = dpx(j,2)                                        ! dpy/dpx1 = dpx/dpy1
                                dpy(j,2) = (px1(j)**2+pz1(j)**2)/gamma(j)**3                      ! dpy/dpy1
                                dpy(j,3) = -1.0D0*(py1(j)*pz1(j))/gamma(j)**3                     ! dpy/dpz1
                                dpz(j,1) = dpx(j,3)                                        ! dpz/dpx1 = dpx/dpz1
                                dpz(j,2) = dpy(j,3)                                        ! dpz/dpy1 = dpy/dpz1
                                dpz(j,3) = (px1(j)**2+py1(j)**2)/gamma(j)**3                      ! dpz/dpz1
                        END DO

                        DO i=1,3
                                dAA(i,k,1,1)=2.0D0*px(1)*a1(k)*dpx(1,i)+2.0D0*qx(1)*b1(k)*dqx(1,i)+&
                                            &2.0D0*rx(1)*c1(k)*drx(1,i)
                                dAA(i,k,1,2)=a1(k)*(dpy(1,i)*px(1)+py(1)*dpx(1,i))+b1(k)*(dqy(1,i)*qx(1)+&
                                            &qy(1)*dqx(1,i))+c1(k)*(dry(1,i)*rx(1)+ry(1)*drx(1,i))
                                dAA(i,k,1,3)=a1(k)*(dpz(1,i)*px(1)+pz(1)*dpx(1,i))+b1(k)*(dqz(1,i)*qx(1)+&
                                            &qz(1)*dqx(1,i))+c1(k)*(drz(1,i)*rx(1)+rz(1)*drx(1,i))
                                dAA(i,k,2,1)=dAA(i,k,1,2)
                                dAA(i,k,2,2)=2.0D0*py(1)*a1(k)*dpy(1,i)+2.0D0*qy(1)*b1(k)*dqy(1,i)+&
                                            &2.0D0*ry(1)*c1(k)*dry(1,i)
                                dAA(i,k,2,3)=a1(k)*(dpy(1,i)*pz(1)+py(1)*dpz(1,i))+b1(k)*(dqy(1,i)*qz(1)+&
                                            &qy(1)*dqz(1,i))+c1(k)*(dry(1,i)*rz(1)+ry(1)*drz(1,i))
                                dAA(i,k,3,1)=dAA(i,k,1,3)
                                dAA(i,k,3,2)=dAA(i,k,2,3)
                                dAA(i,k,3,3)=2.0D0*pz(1)*a1(k)*dpz(1,i)+2.0D0*qz(1)*b1(k)*dqz(1,i)+&
                                            &2.0D0*rz(1)*c1(k)*drz(1,i)
                                dBB(i,k,1,1)=2.0D0*px(2)*a1(k)*dpx(2,i)+2.0D0*qx(2)*b1(k)*dqx(2,i)+&
                                            &2.0D0*rx(2)*c1(k)*drx(2,i)
                                dBB(i,k,1,2)=a1(k)*(dpy(2,i)*px(2)+py(2)*dpx(2,i))+b1(k)*(dqy(2,i)*qx(2)+&
                                            &qy(2)*dqx(2,i))+c1(k)*(dry(2,i)*rx(2)+ry(2)*drx(2,i))
                                dBB(i,k,1,3)=a1(k)*(dpz(2,i)*px(2)+pz(2)*dpx(2,i))+b1(k)*(dqz(2,i)*qx(2)+&
                                            &qz(2)*dqx(2,i))+c1(k)*(drz(2,i)*rx(2)+rz(2)*drx(2,i))
                                dBB(i,k,2,1)=dBB(i,k,1,2)
                                dBB(i,k,2,2)=2.0D0*py(2)*a1(k)*dpy(2,i)+2.0D0*qy(2)*b1(k)*dqy(2,i)+&
                                            &2.0D0*ry(2)*c1(k)*dry(2,i)
                                dBB(i,k,2,3)=a1(k)*(dpy(2,i)*pz(2)+py(2)*dpz(2,i))+b1(k)*(dqy(2,i)*qz(2)+&
                                            &qy(2)*dqz(2,i))+c1(k)*(dry(2,i)*rz(2)+ry(2)*drz(2,i))
                                dBB(i,k,3,1)=dBB(i,k,1,3)
                                dBB(i,k,3,2)=dBB(i,k,2,3)
                                dBB(i,k,3,3)=2.0D0*pz(2)*a1(k)*dpz(2,i)+2.0D0*qz(2)*b1(k)*dqz(2,i)+&
                                            &2.0D0*rz(2)*c1(k)*drz(2,i)
                        END DO

                        DO i=1,3
                                dGu1(i,:,:)=-(1.0D0-lambda(k))*MATMUL(MATMUL(C1inv, dAA(i,k,:,:)),C1inv)
                                dGu2(i,:,:)=-lambda(k)*MATMUL(MATMUL(C1inv, dBB(i,k,:,:)),C1inv)
                        END DO     
            
                        DO j=1,3
                                dXu1(j,:)=MATMUL(dGu1(j,:,:),R1)
                                dFu1(j,:)=MATMUL(R1T,dXu1(j,:))
                                dXu2(j,:)=MATMUL(dGu2(j,:,:),R1)
                                dFu2(j,:)=MATMUL(R1T,dXu2(j,:))
                        END DO
                    
                            DUMMY=lambda(k)*(1.0D0-lambda(k))
                    
!                           dF1(k,7 )=DUMMY*dFu1(1,1)   !dF/dalpha1
!                           dF1(k,8 )=DUMMY*dFu1(2,1)      !dF/dbeta1
!                           dF1(k,9 )=DUMMY*dFu2(1,1)       !dF/dalpha2
!                           dF1(k,10)=DUMMY*dFu2(2,1)   !dF/dbeta2
                    
                        dF1(k,7 )=DUMMY*dFu1(1,1)      !dF/d(px1)
                        dF1(k,8 )=DUMMY*dFu1(2,1)      !dF/d(py1)
                        dF1(k,9 )=DUMMY*dFu1(3,1)      !dF/d(pz1)

                        dF1(k,10)=DUMMY*dFu2(1,1)       !dF/d(px2)
                        dF1(k,11)=DUMMY*dFu2(2,1)       !dF/d(py2)
                        dF1(k,12)=DUMMY*dFu2(3,1)       !dF/d(pz2)
                    
!                        END OF DERIVATIVES of the ellipsoid contact function

         
                            LJ1(k)=sigma0(k)/(r-r/SQRT(F1(k))+sigma0(k))
                            LJ2(k)=LJ1(k)**2
                            LJ3(k)=LJ2(k)*LJ1(k)
                            LJ4(k)=LJ2(k)**2
                            LJ5(k)=LJ4(k)*LJ1(k)
                            LJ6(k)=LJ4(k)*LJ2(k)
                            LJ7(k)=LJ6(k)*LJ1(k)
                            LJ11(k)=LJ5(k)*LJ6(k)
                            LJ12(k)=LJ6(k)**2

! Correction terms to the potential if we require a cutoff at rc. 
                            IF (PARAMONOVCUTOFF) THEN
! r/SQRT(F1(k)) = f(A,B), a parameter ONLY dependent on orientation!

                                       CLJ1(k)=sigma0(k)/(PCUTOFF-r/SQRT(F1(k))+sigma0(k))
                                       CLJ2(k)=CLJ1(k)**2
                                       CLJ3(k)=CLJ2(k)*CLJ1(k)
                                       CLJ4(k)=CLJ2(k)**2
                                       CLJ5(k)=CLJ4(k)*CLJ1(k)
                                       CLJ6(k)=CLJ4(k)*CLJ2(k)
                                       CLJ7(k)=CLJ6(k)*CLJ1(k)
                                       CLJ8(k)=CLJ6(k)*CLJ2(k)
                                       CLJ11(k)=CLJ5(k)*CLJ6(k)
                                       CLJ12(k)=CLJ6(k)**2
                                       CLJ13(k)=CLJ12(k)*CLJ1(k)
                                       CLJ14(k)=CLJ7(k)**2
                            END IF
 
!                        LLJ1=sigma1/(r-r/SQRT(F1(1))+sigma1)
!                        LLJ2=LLJ1**2
!                        LLJ3=LLJ2*LLJ1
!                        LLJ4=LLJ2**2
!                        LLJ5=LLJ4*LLJ1
!                        LLJ6=LLJ4*LLJ2
!                        LLJ7=LLJ6*LLJ1
!                        LLJ11=LLJ5*LLJ6
!                        LLJ12=LLJ6**2
              ! interaction between the extra LJ sites:
                        LLJ1=sigma1/rlj
                        LLJ2=LLJ1**2
                        LLJ3=LLJ2*LLJ1
                        LLJ4=LLJ2**2
                        LLJ5=LLJ4*LLJ1
                        LLJ6=LLJ4*LLJ2
                        LLJ7=LLJ6*LLJ1
                        LLJ11=LLJ5*LLJ6
                        LLJ12=LLJ6**2
 
                            DUMMY=1.0D0/rlj
                            DUMMY=DUMMY**2
                            DO j=1,12
                                dLLJ1(j) =-sigma1*DUMMY*drlj(j)
                            END DO


!                           WRITE(*,*) "dF(lambda(k))=",lambda(k),dF_lambda(k),ITER

                            DUMMY=LJ1(k)/sigma0(k)
                            DUMMY=DUMMY**2
                            DUMMY1=1.0D0/SQRT(F1(k))
                            DUMMY2=DUMMY1**3
 
                            DO j=1,6
                                dLJ1(k,j) =-1.0D0*sigma0(k)*DUMMY*(dr(j)-dr(j)*DUMMY1+r*0.5D0*DUMMY2*dF1(k,j))
                            END DO
        
                            DO j=7,12
                                dLJ1(k,j) =-1.0D0*sigma0(k)*DUMMY*(0.5D0*r*DUMMY2*dF1(k,j))
                            END DO
! inner hard core stuff
!                            DO j=1,6
!                                dLLJ1(j) =-1.0D0*sigma1*DUMMY*(dr(j)-dr(j)*DUMMY1+r*0.5D0*DUMMY2*dF1(1,j))
!                        END DO
        
!                        DO j=7,12
!                                dLLJ1(j) =-1.0D0*sigma1*DUMMY*(0.5D0*r*DUMMY2*dF1(1,j))
!                        END DO


                        IF (PARAMONOVCUTOFF) THEN 
                                DUMMY=CLJ1(k)/sigma0(k)
                                DUMMY=DUMMY**2
                                DUMMY1=1.0D0/SQRT(F1(k))
                                DUMMY2=DUMMY1**3
 
                                DO j=1,6
                                        dCLJ1(k,j) =-1.0D0*sigma0(k)*DUMMY*(-1.0D0*DUMMY1*dr(j)+0.5D0*r*DUMMY2*dF1(k,j))
                                END DO
                                DO j=7,12
                                        dCLJ1(k,j) =-1.0D0*sigma0(k)*DUMMY*(0.5D0*r*DUMMY2*dF1(k,j)) !derivatives wrt to orientation
                                END DO
                        END IF
!                        WRITE(*,*) lambda(k),k,ITER
                  END DO        


!                  WRITE(*,*) 'paramonov.f90> ECF = ', F1(1), F1(2)    
                END IF                         ! if (not(PARAMONOVCUTOFF) or (PARAMONOVCUTOFF and r<PCUTOFF)) then....

        IF (PARAMONOVCUTOFF) THEN
                 IF (r>=PCUTOFF) THEN
                         VDUMMY = 0.0D0
                 ELSE

                  ! work out the spline terms for smooth cutoff of extra LJ sites
                   if(rlj2<ron2) then

                      term2=1.0D0
                      term3=0.0D0
                   else if(rlj2>ron2) then
                          term2=(cut-rlj**2)*(cut-rlj**2)*(cut+2.0D0*rlj**2-3.0D0*ron2)*range2inv3
                          term3=rlj*12.0D0*(cut-rlj**2)*(ron2-rlj**2)*range2inv3 ! d(term2)/dr
                   end if

                   VDUMMY = 0.0D0

!         add corrections to normal LJ pair potential to ensure that potential and 
!         derivatives are zero at rc, and vanish smoothly with no discontinuities. 
!                VDUMMY=epsilon1*(LLJ12-0.0D0*LLJ6) !inner hard core

                VDUMMY=epsilon1*(LLJ12-LLJ6)*term2
!repulsive potential and periodic cutoff corrections
                VDUMMY=VDUMMY+4.0D0*epsilon0*(LJ12(1)+(6.0D0*CLJ12(1))*CLJ2(1)/LJ2(1)-7.0D0*CLJ12(1))
!attractive potential and periodic cutoff corrections
                VDUMMY=VDUMMY+4.0D0*epsilon0*(-LJ6(2)+(-3.0D0*CLJ6(2))*CLJ2(2)/LJ2(2) +4.0D0*CLJ6(2))

                        END IF
        ELSE
!                VDUMMY=(f12 - f6)
!                VDUMMY=epsilon1*(LLJ12-LLJ6)*term2
!                WRITE(*,*) 'VDUMMY for LJ sites',VDUMMY,rlj
                VDUMMY=VDUMMY+4.0D0*epsilon0*(LJ12(1)-LJ6(2))    
        END IF
                      do j=1,12
                        if (PARAMONOVCUTOFF) then
                                if (r>=PCUTOFF) then
                                        dVDUMMY(j)=0.0D0
                                else
                                        dVDUMMY(j)=0.0D0

                                        dVDUMMY(j)=epsilon1*(12.0D0*LLJ11*dLLJ1(j)-6.0D0*LLJ5*dLLJ1(j))*term2 + &
                                                        & epsilon1*(LLJ12-LLJ6)*term3*drlj(j)! extra LJ site derivatives
                                        dVDUMMY(j)=dVDUMMY(j)+4.0D0*epsilon0*(&
                                                &12.0D0*dLJ1(1,j)*(LJ11(1)-CLJ14(1)/LJ3(1))&
                                                &+84.0D0*dCLJ1(1,j)*(CLJ13(1)/LJ2(1)-CLJ11(1))) !repulsive derivative
                                        dVDUMMY(j)=dVDUMMY(j)+4.0D0*epsilon0*(&
                                                &-6.0D0*dLJ1(2,j)*(LJ5(2)-CLJ8(2)/LJ3(2))&
                                                &-24.0D0*dCLJ1(2,j)*(CLJ7(2)/LJ2(2)-CLJ5(2))) !attractive derivative


!4.0D0*epsilon0*(-6.0D0*dLJ1(2,j)*LJ5(2)+24.0D0*dCLJ1(2,j)*CLJ5(2)-24.0D0*dCLJ1(2,j)*CLJ7(2)/LJ2(2)+6.0D0*dLJ1(2,j)*CLJ8(2)/LJ3(2)) 
!-CLJ8(2)/LJ3(2))-24.0D0*dCLJ1(2,j)*(CLJ7(2)/LJ2(2)))

                                end if
                        else
!                                 dVDUMMY(j)=epsilon1*(12.0D0*LLJ11*dLLJ1(j)-6.0D0*LLJ5*dLLJ1(j))&
!                                              &+4.0D0*epsilon0*(12.0D0*LJ11(1)*dLJ1(1,j)-6.0D0*LJ5(2)*dLJ1(2,j))
                                  dVDUMMY(j)=4.0D0*epsilon0*(12.0D0*LJ11(1)*dLJ1(1,j)-6.0D0*LJ5(2)*dLJ1(2,j))

                        end if
!                Write(*,*) j, dVDUMMY(j)
                end do !j=1,12
!       WRITE(*,*) 'i,j,dVDUMMY', J1,J3,alpha(2),da(2,:)
!       write(*,*) pz1(2),dVDUMMY(1:12)
!                END DO
!                WRITE(*,*) dCLJ1(1,:),dCLJ1(2,:)

!                        z1 = z1+0.0000001D0
!                        Fcheck(l)=VDUMMY
!                        IF(l==2) Dcheck=dVDUMMY(3)
!                        IF(l==3) THEN
!                         WRITE(*,*) (Fcheck(3)-Fcheck(1))/0.0000002D0, Dcheck
!                         STOP
!                        END IF
!             END DO    ! for derivative check


                EGB=EGB+VDUMMY
    
                VT(J1) = VT(J1) + VDUMMY
                VT(J3) = VT(J3) + VDUMMY        ! pair potentials

                V(J2-2)=V(J2-2)+dVDUMMY(1)
                V(K2-2)=V(K2-2)+dVDUMMY(4)
                V(J2-1)=V(J2-1)+dVDUMMY(2)
                V(K2-1)=V(K2-1)+dVDUMMY(5)
                V(J2  )=V(J2  )+dVDUMMY(3)        
                V(K2  )=V(K2  )+dVDUMMY(6)       
                V(3*REALNATOMS+J2-2)=V(3*REALNATOMS+J2-2)+dVDUMMY(7)        
                V(3*REALNATOMS+K2-2)=V(3*REALNATOMS+K2-2)+dVDUMMY(10)        
                V(3*REALNATOMS+J2-1)=V(3*REALNATOMS+J2-1)+dVDUMMY(8)        
                V(3*REALNATOMS+K2-1)=V(3*REALNATOMS+K2-1)+dVDUMMY(11)       
                V(3*REALNATOMS+J2  )=V(3*REALNATOMS+J2  )+dVDUMMY(9)        
                V(3*REALNATOMS+K2  )=V(3*REALNATOMS+K2  )+dVDUMMY(12)       

!                WRITE(*,*) 'Dumping derivatives:',J1,J3
!                WRITE(*,*) dVDUMMY(:)
!                WRITE(*,*) X(:)
!                =SQRT((x1-x2)**2,(y1-y2)**2,(z1-z2)**2)

! some debug printing

!WRITE(*,*) 'x1,y1,z1'
!WRITE(*,*) x1,y1,z1
!WRITE(*,*) 'x2,y2,z2'
!WRITE(*,*) x2,y2,z2
!WRITE(*,*) 'px1,py1,pz1'
!WRITE(*,*) px1(1),py1(1),pz1(1)
!WRITE(*,*) 'px2,py2,pz2'
!WRITE(*,*) px1(2),py1(2),pz1(2)

!WRITE(*,*) 'energy:' , VDUMMY
!WRITE(*,*) 'lambda, repulsive contact function', lambda(1), F1(1) 
!WRITE(*,*) 'lambda, attractive contact function', lambda(2), F1(2)
!STOP
            END DO      !DO J3=J1+1,REALNATOMS
           END DO       !DO J1=1,REALNATOMS
        
!           IF(STEST) CALL PARAMONOVSECDER(X,.TRUE.)

!      WRITE(*,*) "Exiting energy cycle"
999    CONTINUE      
       END SUBROUTINE OLDPARAMONOV


SUBROUTINE PARAMONOVNUMFIRSTDER(OLDX,STEST)
use commons, only : natoms
!use key
use modhess
implicit none

DOUBLE PRECISION  :: V(3*NATOMS),ENERGY,X(3*NATOMS),NUMGRAD(3*NATOMS),TRUEGRAD(3*NATOMS),&
        &            OLDX(3*NATOMS),ETEMP(2,3*NATOMS),ksi
INTEGER    :: i,j
LOGICAL    :: GTEST,STEST


ksi=0.0001D0
X(:)=OLDX(:)

!         CALL OLDPARAMONOV(X,TRUEGRAD,ENERGY,GTEST,.FALSE.)
         CALL PYGPERIODIC(X,TRUEGRAD,ENERGY,.TRUE.)


!WRITE(*,*) 'sf344> in Paramonovsecder'
DO i=1,3*NATOMS

         X(i)=X(i)-ksi

!         CALL OLDPARAMONOV(X,V,ENERGY,GTEST,.FALSE.)
         CALL PYGPERIODIC(X,V,ENERGY,.TRUE.)

         ETEMP(1,i)=ENERGY

         X(i)=X(i)+2.0D0*ksi

!         CALL OLDPARAMONOV(X,V,ENERGY,GTEST,.FALSE.)
         CALL PYGPERIODIC(X,V,ENERGY,.TRUE.)
    
         ETEMP(2,i)=ENERGY
         NUMGRAD(i)=(ETEMP(2,i)-ETEMP(1,i))/(2.0D0*ksi)

         X(i)=X(i)-ksi

WRITE(*,'(3F20.10)') TRUEGRAD(i),NUMGRAD(i), 100.0D0*(TRUEGRAD(i)-NUMGRAD(i))/TRUEGRAD(i)

END DO

STOP
!WRITE(*,*) 'sf344> exiting AMBERsecder'
!WRITE(*,*) 'HESSIAN:'
!WRITE(*,'12F10.3') HESS(:,:)
END SUBROUTINE PARAMONOVNUMFIRSTDER





SUBROUTINE ECF(OVERLAPTEST,ECFvalue,x1,x2,y1,y2,z1,z2,px11,px12,py11,py12,pz11,pz12,a,b,c)


!Calculates the value of the ellipsoid contact function

!use commons
IMPLICIT NONE

DOUBLE PRECISION, intent(in)       :: x1,x2,y1,y2,z1,z2,px11,px12,py11,py12,pz11,pz12,a,b,c
DOUBLE PRECISION, intent(out)      :: ECFvalue
DOUBLE PRECISION                   :: a1(2),b1(2),c1(2),lambdaold,dF_lambdaold,&
                &       DUMMY,DUMMY1,DUMMY2,S_DUMMY1,C_DUMMY1,S_DUMMY2,C_DUMMY2,px1(2),py1(2),pz1(2),&
                &       alpha(2),beta(2),gamma(2),qx(2),qy(2),qz(2),rx(2),ry(2),rz(2),px(2),py(2),pz(2),&
                &       sa(2),ca(2),sb(2),cb(2),sg(2),cg(2)
DOUBLE PRECISION  :: S_lambda(1),F1(1),lambda(2), A11(3,3),B11(3,3),C_1(6),C1inv(3,3),W(3),R1(3),R1T(1,3),&
                &       dF_lambda(1),G_invR(3),G_invRT(1,3),DPRAND, RANDOM, STEPSIZE, CONVG, r, PI
INTEGER    :: i,j,k,ITER
!LOGICAL, intent(out)    :: OVERLAPTEST
LOGICAL   :: OVERLAPTEST

a1(1) = a**2
b1(1) = b**2
c1(1) = c**2

px1(1) = px11
px1(2) = px12
py1(1) = py11
py1(2) = py12
pz1(1) = pz11
pz1(2) = pz12

PI = ATAN(1.0)*4.0D0

k=1
        r=SQRT((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
                   
 R1(1)=x2-x1
 R1(2)=y2-y1
 R1(3)=z2-z1
 
 R1T(1,1)=R1(1)
 R1T(1,2)=R1(2)
 R1T(1,3)=R1(3)

           stepsize=1.0D0
            lambda(1)=0.5D0
            ITER=0
            convg=1.0D-6
            dF_lambdaold=1.0D6
 
   DO i=1,2
     gamma(i) = sqrt(px1(i)**2+py1(i)**2+pz1(i)**2)

     px(i) = px1(i)/gamma(i)
     py(i) = py1(i)/gamma(i)
     pz(i) = pz1(i)/gamma(i)

     IF(px(i).eq.0.0d0) THEN
       IF(py(i)>=0.0D0) THEN
        alpha(i) = PI/2           ! Euler angle alpha 
       ELSE
        alpha(i) = -1.0D0*PI/2
       END IF
     ELSE
        IF(py(i)>=0.0D0) THEN
          IF(px(i)>0.0D0) THEN
            alpha(i) = 1.0D0*atan(py(i)/px(i))   ! first quadrant
          ELSE    ! px<0
            alpha(i) = 1.0D0*atan(py(i)/px(i))+PI       ! should be in the second quadrant
          END IF
        ELSE IF(py(i)<0.0D0) THEN
          IF(px(i)>0.0D0) THEN
            alpha(i) = 1.0D0*atan(py(i)/px(i))    ! fourth quadrant
          ELSE    ! px<0
            alpha(i) = 1.0D0*atan(py(i)/px(i))-PI            ! third quadrant
          END IF
        END IF 
     END IF
 
        beta(i) = 1.0D0*acos(pz(i))

        DUMMY1 = pz1(i)/(sqrt(px1(i)**2+py1(i)**2)*gamma(i)**2)
        
     sa(i) = sin(alpha(i))
     sb(i) = sin(beta(i))
     sg(i) = sin(gamma(i))
     ca(i) = cos(alpha(i))
     cb(i) = cos(beta(i))
     cg(i) = cos(gamma(i))

     qx(i) = -sg(i)*cb(i)*ca(i)-cg(i)*sa(i)
     qy(i) = cg(i)*ca(i)-sg(i)*cb(i)*sa(i)
     qz(i) = sg(i)*sb(i)
     rx(i) = cg(i)*cb(i)*ca(i)-sg(i)*sa(i)
     ry(i) = cg(i)*cb(i)*sa(i)+sg(i)*ca(i)
     rz(i) = -cg(i)*sb(i)
   END DO  ! i=1,2

      A11(1,1)=px(1)**2*a1(k)+qx(1)**2*b1(k)+rx(1)**2*c1(k)
      A11(1,2)=py(1)*px(1)*a1(k)+qy(1)*qx(1)*b1(k)+ry(1)*rx(1)*c1(k)
      A11(1,3)=px(1)*pz(1)*a1(k)+qx(1)*qz(1)*b1(k)+rx(1)*rz(1)*c1(k)
      A11(2,1)=A11(1,2)
      A11(2,2)=py(1)**2*a1(k)+qy(1)**2*b1(k)+ry(1)**2*c1(k)
      A11(2,3)=py(1)*pz(1)*a1(k)+qy(1)*qz(1)*b1(k)+ry(1)*rz(1)*c1(k)
      A11(3,1)=A11(1,3)
      A11(3,2)=A11(2,3)
      A11(3,3)=pz(1)**2*a1(k)+qz(1)**2*b1(k)+rz(1)**2*c1(k)

      B11(1,1)=px(2)**2*a1(k)+qx(2)**2*b1(k)+rx(2)**2*c1(k)
      B11(1,2)=py(2)*px(2)*a1(k)+qy(2)*qx(2)*b1(k)+ry(2)*rx(2)*c1(k)
      B11(1,3)=px(2)*pz(2)*a1(k)+qx(2)*qz(2)*b1(k)+rx(2)*rz(2)*c1(k)
      B11(2,1)=B11(1,2)
      B11(2,2)=py(2)**2*a1(k)+qy(2)**2*b1(k)+ry(2)**2*c1(k)
      B11(2,3)=py(2)*pz(2)*a1(k)+qy(2)*qz(2)*b1(k)+ry(2)*rz(2)*c1(k)
      B11(3,1)=B11(1,3)
      B11(3,2)=B11(2,3)
      B11(3,3)=pz(2)**2*a1(k)+qz(2)**2*b1(k)+rz(2)**2*c1(k)


           DO
                  RANDOM=DPRAND()
                   dF_lambda(:)=0.0D0
                   W(:)=0.0D0

                  C_1(1)=(1.0D0-lambda(k))*A11(1,1)+lambda(k)*B11(1,1)
                  C_1(2)=(1.0D0-lambda(k))*A11(1,2)+lambda(k)*B11(1,2)
                  C_1(3)=(1.0D0-lambda(k))*A11(1,3)+lambda(k)*B11(1,3)
                  C_1(4)=(1.0D0-lambda(k))*A11(2,2)+lambda(k)*B11(2,2)
                  C_1(5)=(1.0D0-lambda(k))*A11(2,3)+lambda(k)*B11(2,3)
                  C_1(6)=(1.0D0-lambda(k))*A11(3,3)+lambda(k)*B11(3,3)

 
 
        ! WRITE(*,*) 'Calling SVERT, matrix to be inverted:'
        ! WRITE(*,*) C_1(:)
        ! WRITE(*,*) ITER
        ! WRITE(*,*) W(:)
          
             CALL SVERT(C_1,3,W)  !inverting the symmetric matrix C1 (or G, eq.2.6 in ECF paper)
     

      
             C1inv(1,1)=C_1(1)
             C1inv(1,2)=C_1(2)
             C1inv(1,3)=C_1(3)
             C1inv(2,1)=C1inv(1,2)
             C1inv(2,2)=C_1(4)
             C1inv(2,3)=C_1(5)
             C1inv(3,1)=C1inv(1,3)
             C1inv(3,2)=C1inv(2,3)
             C1inv(3,3)=C_1(6)


        ! WRITE(*,*) 'Inverted matrix:'
        ! WRITE(*,*) C1inv(:,:)
     
            G_invR=MATMUL(C1inv,R1)
     
             G_invRT(1,1)=G_invR(1)
             G_invRT(1,2)=G_invR(2)
             G_invRT(1,3)=G_invR(3)
     
    
      !derivative of the contact function with respect to lambda
       dF_lambda=MATMUL(G_invRT,MATMUL((1.0D0-lambda(k))**2*A11-lambda(k)**2*B11,G_invR))   
          DUMMY=dF_lambda(1)
          DUMMY1=lambda(k)

             IF(ABS(dF_lambda(1))<convg) THEN
                EXIT
             ELSE
                ITER=ITER+1

        ! if the new derivative is larger than the old one, reset lambda to its previous value and decrease stepsize
                IF(ABS(dF_lambda(1))>ABS(dF_lambdaold)) THEN
!                WRITE(*,*) 'dF_lambda>df_lambdaold', dF_lambda, df_lambdaold
                    lambda(k)=lambdaold
                    dF_lambda(1)=dF_lambdaold
                    !IF(stepsize>=0.000000001D0)
                    stepsize=stepsize/3.0D0
                END IF

                IF(lambda(k)<0.OR.lambda(k)>1) THEN
                   WRITE(*,*) "lambda value not in [0,1], resetting"
                    stepsize=1.0D0
                    lambda(k)=0.5
                    dF_lambda(1)=0.1D0
                    DUMMY=dF_lambda(1)
                    DUMMY1=lambda(k)
                END IF

                IF(ABS(dF_lambda(1))<=0.1D0) THEN
                    lambda(k)=lambda(k)+dF_lambda(1)*stepsize*RANDOM
                ELSE IF(ABS(dF_lambda(1))>1.0D0) THEN
                    lambda(k)=lambda(k)+(dF_lambda(1)/ABS(dF_lambda(1)))*RANDOM*stepsize*0.2
                ELSE
                    lambda(k)=lambda(k)+dF_lambda(1)*stepsize*0.1*RANDOM
                END IF
 
  lambdaold=DUMMY1
  dF_lambdaold=DUMMY

             END IF
!           IF(ITER.GT.1000) THEN
!
!                WRITE(*,*) lambda(k),dF_lambda(1),ITER,px,py,py/px,alpha
!                WRITE(*,*) x1,y1,z1
!                WRITE(*,*) x2,y2,z2
!                WRITE(*,*) px1,py1,pz1
!                WRITE(*,*) px2,py2,pz2
!                WRITE(*,*) SQRT(px1**2+py1**2+pz1**2)
!                WRITE(*,*) SQRT(px2**2+py2**2+pz2**2)
!           END IF
           END DO
   
           S_lambda=lambda(k)*(1.0D0-lambda(k))*MATMUL(R1T,G_invR)  ! this is the ellipsoid contact function
           F1(k)=S_lambda(1)
ECFvalue=S_lambda(1)

IF(ECFvalue<1.0D0) THEN
 OVERLAPTEST=.TRUE.
ELSE
 OVERLAPTEST=.FALSE.
END IF

!WRITE(*,*) k
!WRITE(*,*) OVERLAPTEST

END SUBROUTINE ECF


SUBROUTINE OLDECF2(OVERLAPTEST,ECFvalue,x1,x2,y1,y2,z1,z2,px1,px2,py1,py2,pz1,pz2,a,b,c)


!Calculates the value of the ellipsoid contact function

!use commons
IMPLICIT NONE

DOUBLE PRECISION, intent(in)       :: x1,x2,y1,y2,z1,z2,px1,px2,py1,py2,pz1,pz2,a,b,c
DOUBLE PRECISION, intent(out)      :: ECFvalue
DOUBLE PRECISION                   :: a1(2),b1(2),c1(2),lambdaold,dF_lambdaold,&
                &       DUMMY,DUMMY1,DUMMY2,S_DUMMY1,C_DUMMY1,S_DUMMY2,C_DUMMY2
DOUBLE PRECISION  :: S_lambda(1),F1(1),lambda(2), A11(3,3),B11(3,3),C_1(6),C1inv(3,3),W(3),R1(3),R1T(1,3),&
                &       dF_lambda(1),G_invR(3),G_invRT(1,3),DPRAND, RANDOM, STEPSIZE, CONVG, r, PI
INTEGER    :: i,j,k,ITER
LOGICAL, intent(out)    :: OVERLAPTEST

a1(1) = a**2
b1(1) = b**2
c1(1) = c**2

PI = ATAN(1.0)*4.0D0

k=1
        r=SQRT((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
                   
 R1(1)=x2-x1
 R1(2)=y2-y1
 R1(3)=z2-z1
 
 R1T(1,1)=R1(1)
 R1T(1,2)=R1(2)
 R1T(1,3)=R1(3)

!       dr(1)=-1.0D0/r*(x2-x1)
!       dr(2)=-1.0D0/r*(y2-y1)
!       dr(3)=-1.0D0/r*(z2-z1)
!       dr(4)=-dr(1)
!       dr(5)=-dr(2)
!       dr(6)=-dr(3)
                          
!      px1=C_alpha1*C_beta1
!      py1=C_alpha1*S_beta1 
!      pz1=S_alpha1        
!      px2=C_alpha2*C_beta2
!      py2=C_alpha2*S_beta2
!      pz2=S_alpha2

           stepsize=1.0D0
            lambda(k)=0.5D0
            ITER=0
            convg=1.0D-2
            dF_lambdaold=1.0D6
 

     IF(px1.eq.0.0d0) THEN
        DUMMY1 = PI/2
     ELSE
        DUMMY1 = atan(py1/px1)
     END IF

     IF(px2.eq.0.0d0) THEN
        DUMMY2 = PI/2
     ELSE
        DUMMY2 = atan(py2/px2)
     END IF

        S_DUMMY1 = sin(DUMMY1)
        C_DUMMY1 = cos(DUMMY1)
        S_DUMMY2 = sin(DUMMY2)
        C_DUMMY2 = cos(DUMMY2)

   
      DUMMY = SQRT(px1**2+py1**2+pz1**2)

      A11(1,1)=(px1/DUMMY)**2*a1(k)+S_DUMMY1**2*b1(k)+C_DUMMY1**2*(pz1/DUMMY)**2*c1(k)
      A11(1,2)=py1*px1*a1(k)/(DUMMY**2)-C_DUMMY1*S_DUMMY1*b1(k)+pz1**2*C_DUMMY1*S_DUMMY1*c1(k)/(DUMMY**2)
      A11(1,3)=px1*pz1*(a1(k)-c1(k))/(DUMMY**2)
      A11(2,1)=A11(1,2)
      A11(2,2)=(py1/DUMMY)**2*a1(k)+C_DUMMY1**2*b1(k)+(pz1/DUMMY)**2*S_DUMMY1**2*c1(k)
      A11(2,3)=py1*pz1*(a1(k)-c1(k))/(DUMMY**2)
      A11(3,1)=A11(1,3)
      A11(3,2)=A11(2,3)
      A11(3,3)=(pz1/DUMMY)**2*a1(k)+(px1**2+py1**2)*c1(k)/(DUMMY**2)

      DUMMY = SQRT(px2**2+py2**2+pz2**2)

      B11(1,1)=(px2/DUMMY)**2*a1(k)+S_DUMMY2**2*b1(k)+C_DUMMY2**2*(pz2/DUMMY)**2*c1(k)
      B11(1,2)=py2*px2*a1(k)/(DUMMY**2)-C_DUMMY2*S_DUMMY2*b1(k)+pz2**2*C_DUMMY2*S_DUMMY2*c1(k)/(DUMMY**2)
      B11(1,3)=px2*pz2*(a1(k)-c1(k))/(DUMMY**2)
      B11(2,1)=B11(1,2)
      B11(2,2)=(py2/DUMMY)**2*a1(k)+C_DUMMY2**2*b1(k)+(pz2/DUMMY)**2*S_DUMMY2**2*c1(k)
      B11(2,3)=py2*pz2*(a1(k)-c1(k))/(DUMMY**2)
      B11(3,1)=B11(1,3)
      B11(3,2)=B11(2,3)
      B11(3,3)=(pz2/DUMMY)**2*a1(k)+(px2**2+py2**2)*c1(k)/(DUMMY**2)


           DO
                  RANDOM=DPRAND()
                   dF_lambda(:)=0.0D0
                   W(:)=0.0D0

                  C_1(1)=(1.0D0-lambda(k))*A11(1,1)+lambda(k)*B11(1,1)
                  C_1(2)=(1.0D0-lambda(k))*A11(1,2)+lambda(k)*B11(1,2)
                  C_1(3)=(1.0D0-lambda(k))*A11(1,3)+lambda(k)*B11(1,3)
                  C_1(4)=(1.0D0-lambda(k))*A11(2,2)+lambda(k)*B11(2,2)
                  C_1(5)=(1.0D0-lambda(k))*A11(2,3)+lambda(k)*B11(2,3)
                  C_1(6)=(1.0D0-lambda(k))*A11(3,3)+lambda(k)*B11(3,3)

 
 
        ! WRITE(*,*) 'Calling SVERT, matrix to be inverted:'
        ! WRITE(*,*) C_1(:)
        ! WRITE(*,*) ITER
        ! WRITE(*,*) W(:)
          
             CALL SVERT(C_1,3,W)  !inverting the symmetric matrix C1 (or G, eq.2.6 in ECF paper)
     

      
             C1inv(1,1)=C_1(1)
             C1inv(1,2)=C_1(2)
             C1inv(1,3)=C_1(3)
             C1inv(2,1)=C1inv(1,2)
             C1inv(2,2)=C_1(4)
             C1inv(2,3)=C_1(5)
             C1inv(3,1)=C1inv(1,3)
             C1inv(3,2)=C1inv(2,3)
             C1inv(3,3)=C_1(6)


        ! WRITE(*,*) 'Inverted matrix:'
        ! WRITE(*,*) C1inv(:,:)
     
            G_invR=MATMUL(C1inv,R1)
     
             G_invRT(1,1)=G_invR(1)
             G_invRT(1,2)=G_invR(2)
             G_invRT(1,3)=G_invR(3)
     
    
     
             dF_lambda=MATMUL(G_invRT,MATMUL((1.0D0-lambda(k))**2*A11-lambda(k)**2*B11,G_invR)) !derivative of the contact function with respect to lambda
             DUMMY=dF_lambda(1)
             DUMMY1=lambda(k)
! WRITE(*,*) "Starting subroutine Paramonov",ITER, lambda(k), dF_lambda(:)
! WRITE(*,*) G_invRT(1,:),MATMUL((1.0D0-lambda(k))**2*A11-lambda(k)**2*B11,G_invR)

             IF(ABS(dF_lambda(1))<convg) THEN
                  EXIT
             ELSE
                 ITER=ITER+1

        ! if the new derivative is larger than the old one, reset lambda to its previous value and decrease stepsize
                IF(ABS(dF_lambda(1))>ABS(dF_lambdaold)) THEN
                    lambda(k)=lambdaold
                    dF_lambda(1)=dF_lambdaold
                    !IF(stepsize>=0.000000001D0)
                    stepsize=stepsize/3.0D0
                END IF

                IF(lambda(k)<0.OR.lambda(k)>1) THEN
!                   WRITE(*,*) "lambda value not in [0,1], resetting"
                    stepsize=1.0D0
                    lambda(k)=0.5
                    dF_lambda(1)=0.1D0
                    DUMMY=dF_lambda(1)
                    DUMMY1=lambda(k)
                END IF


                IF(ABS(dF_lambda(1))<=0.1D0) THEN
                    lambda(k)=lambda(k)+dF_lambda(1)*stepsize*RANDOM
                ELSE IF(ABS(dF_lambda(1))>1.0D0) THEN
                    lambda(k)=lambda(k)+(dF_lambda(1)/ABS(dF_lambda(1)))*RANDOM*stepsize*0.2
                ELSE
                    lambda(k)=lambda(k)+dF_lambda(1)*stepsize*0.1*RANDOM
                END IF

 
  lambdaold=DUMMY1
  dF_lambdaold=DUMMY

  
             END IF
!           IF(ITER.GT.1000) THEN

!                WRITE(*,*) lambda(k),dF_lambda(1),ITER
!                WRITE(*,*) x1,y1,z1
!                WRITE(*,*) x2,y2,z2
!                WRITE(*,*) px1,py1,pz1
!                WRITE(*,*) px2,py2,pz2
!                WRITE(*,*) SQRT(px1**2+py1**2+pz1**2)
!                WRITE(*,*) SQRT(px2**2+py2**2+pz2**2)
!           END IF
           END DO
   
    
            S_lambda=lambda(k)*(1.0D0-lambda(k))*MATMUL(R1T,G_invR)  ! this is the ellipsoid contact function
            F1(k)=S_lambda(1)

ECFvalue=S_lambda(1)

IF(ECFvalue<1.0D0) THEN
 OVERLAPTEST=.TRUE.
ELSE
 OVERLAPTEST=.FALSE.
END IF

!WRITE(*,*) k
!WRITE(*,*) OVERLAPTEST

END SUBROUTINE OLDECF2



SUBROUTINE GAYBERNE(X,V,EGB,GTEST,STEST)
use commons
implicit none

! Calculation of the Gay-Berne potential and derivatives
!
! Variables: x1, y1, z1, x2, y2, z2, px1, py1, pz1, px2, py2, pz2 
! V: 1D array holding the derivatives

LOGICAL   :: GTEST, STEST, OVERLAP, OVERLAPTEST
LOGICAL   :: OVERLAPTESTT
INTEGER   :: i,J1,J2,J3,K1, K2, REALNATOMS
DOUBLE PRECISION :: X(3*NATOMS), V(3*NATOMS), VD(3*NATOMS), x1, y1, z1, x2, y2, z2, px1, py1, pz1, px2, py2, pz2,&
&    ru1, ru2, r, r2, u1u2, EGB, DPRAND,DERIV, NUMDER(4),DUMMY11, DUMMY22, DUMMY33,&
&    khi, khi1, sigma, sigmacut, FLJ, FLJc, FLJshift, epsilonnu, epsilonmu, epsilon, DUMMY, DUMMY1, &
&    DUMMY2, DUMMY3,DUMMY4, DUMMY5, VDUMMY, RANDOM, &
&    LJ1, LJ2, LJ3, LJ4, LJ6, LJ7, LJ12, LJ13, ECFvalue, &
&    LJ1c, LJ2c, LJ3c, LJ4c, LJ6c, LJ7c, LJ12c, LJ13c, &
&    LJ1cc, LJ2cc, LJ3cc, LJ4cc, LJ6cc, LJ7cc, LJ12cc, LJ13cc, &
&    rcut, rcutu1, rcutu2, Vcut, &
&    dr_x1, dr_x2, dr_y1, dr_y2, dr_z1, dr_z2, &
&    dV_px1, dFLJ_px1,  dFLJc_px1, dsigma_px1, depsilonnu_px1, depsilonmu_px1, depsilon_px1, dru1_px1, du1u2_px1, &
&    dV_px2, dFLJ_px2,  dFLJc_px2, dsigma_px2, depsilonnu_px2, depsilonmu_px2, depsilon_px2, dru2_px2, du1u2_px2, &
&    dV_py1, dFLJ_py1,  dFLJc_py1, dsigma_py1, depsilonnu_py1, depsilonmu_py1, depsilon_py1, dru1_py1, du1u2_py1, &
&    dV_py2, dFLJ_py2,  dFLJc_py2, dsigma_py2, depsilonnu_py2, depsilonmu_py2, depsilon_py2, dru2_py2, du1u2_py2, &
&    dV_pz1, dFLJ_pz1,  dFLJc_pz1, dsigma_pz1, depsilonnu_pz1, depsilonmu_pz1, depsilon_pz1, dru1_pz1, du1u2_pz1, &
&    dV_pz2, dFLJ_pz2,  dFLJc_pz2, dsigma_pz2, depsilonnu_pz2, depsilonmu_pz2, depsilon_pz2, dru2_pz2, du1u2_pz2, &
&    dV_x1 , dFLJ_x1 ,  dFLJc_x1 , dsigma_x1, depsilonnu_x1, depsilonmu_x1, depsilon_x1, dru1_x1, dru2_x1,  &
&    dV_x2 , dFLJ_x2 ,  dFLJc_x2 , dsigma_x2, depsilonnu_x2, depsilonmu_x2, depsilon_x2, dru1_x2, dru2_x2,  &
&    dV_y1 , dFLJ_y1 ,  dFLJc_y1 , dsigma_y1, depsilonnu_y1, depsilonmu_y1, depsilon_y1, dru1_y1, dru2_y1,  &
&    dV_y2 , dFLJ_y2 ,  dFLJc_y2 , dsigma_y2, depsilonnu_y2, depsilonmu_y2, depsilon_y2, dru1_y2, dru2_y2,  &
&    dV_z1 , dFLJ_z1 ,  dFLJc_z1 , dsigma_z1, depsilonnu_z1, depsilonmu_z1, depsilon_z1, dru1_z1, dru2_z1,  &
&    dV_z2 , dFLJ_z2 ,  dFLJc_z2 , dsigma_z2, depsilonnu_z2, depsilonmu_z2, depsilon_z2, dru1_z2, dru2_z2,  &
&    kappa,kappa1     ! sigmaPAR/sigmaPERP, Es/Ee

DOUBLE PRECISION :: sigma0, epsilon0, mu, nu,  &  ! Default parametres of the Gay-Berne potential:
&    sigmaPAR, sigmaPERP,Ee,Es                                ! SigmaPAR=3.0, mu=2.0, nu=1.0
                                                              ! Es=1.0, Ee=0.2: strength parameters for side-by-side
                                                              ! and end-to-end configurations

! WRITE(*,*) "GB parametres:", GBANISOTROPYR, GBWELLDEPTHR, PSIGMA0, PEPSILON0
! WRITE(*,*) 'Gay-Berne routine called', OVERLAPTESTT

IF (OVERLAPTESTT) THEN
    OVERLAPTEST = .TRUE.
ELSE
    OVERLAPTEST = .FALSE.
END IF
sigma0=PSIGMA0(1)
epsilon0=PEPSILON0
mu = GBMU
nu = GBNU


sigmaPAR=GBANISOTROPYR
sigmaPERP=1.0D0

Ee=GBWELLDEPTHR
Es=1.0D0

!a=GBANISOTROPYR/2.0D0
!b=0.5D0

!WRITE(*,*) Ee, Es
!WRITE(*,*) sigmaPAR,sigmaPERP

! initialization of parametres

    khi=(sigmaPAR**2-sigmaPERP**2)/(sigmaPAR**2+sigmaPERP**2)    ! anisotropy parameter of ellipsoid
    khi1=(Es**(1.0D0/mu)-Ee**(1.0D0/mu))/(Es**(1.0D0/mu)+Ee**(1.0D0/mu))   ! another anisotropy-like parameter

    kappa=sigmaPAR/sigmaPERP
    kappa1=Es/Ee

    rcut = (GBANISOTROPYR+2.0D0)*sigma0!*2.0D0
REALNATOMS=NATOMS/2
i=REALNATOMS

111 CONTINUE

!VT(1:REALNATOMS)=0.0D0
V(1:6*REALNATOMS)=0.0D0
EGB=0.0D0

!WRITE(*,*) 'REALNATOMS', REALNATOMS




DO J1=1,REALNATOMS

    J2=3*J1

    x1 = X(J2-2)
    y1 = X(J2-1)
    z1 = X(J2  )

 ! px, py, pz: angle-axis coordinate system, trying to convert the potential to these ones

    px1 = X(3*REALNATOMS+J2-2)
    py1 = X(3*REALNATOMS+J2-1)
    pz1 = X(3*REALNATOMS+J2  )

!try to normalize px, py, pz

        DUMMY = SQRT(px1**2+py1**2+pz1**2)
        px1 = px1 / DUMMY
        py1 = py1 / DUMMY
        pz1 = pz1 / DUMMY

! first make pz1 always positive (doesn't matter for uniaxial ellipsoids, vec1 = -vec1)

!    IF (pz1.LT.0.0D0) THEN
!        px1 = -px1
!        py1 = -py1
!        pz1 = -pz1
!    END IF


! put the normalized coordinates back into the coordinate array
    X(3*REALNATOMS+J2-2) = px1
    X(3*REALNATOMS+J2-1) = py1
    X(3*REALNATOMS+J2  ) = pz1


    DO J3=J1+1,REALNATOMS
  
      VDUMMY=0.0D0
         K2=3*J3
         x2 = X(K2-2)
         y2 = X(K2-1)
         z2 = X(K2  )
         px2 = X(3*REALNATOMS+K2-2)
         py2 = X(3*REALNATOMS+K2-1)
         pz2 = X(3*REALNATOMS+K2  )

!! inner loop to check numerical derivatives

!DO i=0,3

!py2= py2+0.0001D0
!IF (i==3) py2 = py2-0.0004D0


!try to normalize px, py, pz

!IF(i==0) THEN
        DUMMY = SQRT(px2**2+py2**2+pz2**2)
        px2 = px2 / DUMMY
        py2 = py2 / DUMMY
        pz2 = pz2 / DUMMY
!END IF


! first make pz2 always positive (doesn't matter for uniaxial ellipsoids, vec1 = -vec1)

!    IF (pz2.LT.0.0D0) THEN
!        px2 = -px2
!        py2 = -py2
!        pz2 = -pz2
!    END IF

! put the normalized coordinates back into the coordinate array
         X(3*REALNATOMS+K2-2) = px2
         X(3*REALNATOMS+K2-1) = py2
         X(3*REALNATOMS+K2  ) = pz2


!         K2=3*k
!         x2=X(K2-2)
!         y2=X(K2-1)
!         z2=X(K2  )
!         alpha2=X(3*REALNATOMS+K2-2)
!         beta2=X(3*REALNATOMS+K2-1)
        
! alphaR=ATAN((z2-z1)/SQRT((x2-x1)**2+(y2-y1)**2))          !alphaR: angle of the r vector and xy plane

        r=SQRT((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)

! use cutoff
!IF (r>rcut) then
!         VDUMMY = 0.0D0

!         dV_x1=0.0D0
!         dV_x2=0.0D0
!         dV_y1=0.0D0
!         dV_y2=0.0D0
!         dV_z1=0.0D0
!         dV_z2=0.0D0

!         dV_px1=0.0D0
!         dV_py1=0.0D0
!         dV_pz1=0.0D0

!         dV_px2=0.0D0
!         dV_py2=0.0D0
!         dV_pz2=0.0D0
!WRITE(*,*) r, rcut
!STOP

!ELSE

! check for cold fusion
!      IF(ECFVALUE<0.2D0)      THEN   
!         WRITE(*,*) 'sf344> cold fusion detected', ECFvalue 
!         VDUMMY=-1.0D20
!         EGB=-1.0D20
!         V(:)=0.0D0
!         GOTO 999
!      END IF

     
!        IF(OVERLAPTEST) THEN
!                WRITE(*,*) 'gay-berne.f90> overlap detected', ECFvalue
!        END IF

!      WRITE(*,*) 'Value of the ECF:', ECFvalue, OVERLAPTEST

!        WRITE(*,*) 'Checking values:', px2,py2,pz2,i


        dr_x1 = (x1-x2)/r
        dr_y1 = (y1-y2)/r
        dr_z1 = (z1-z2)/r
        dr_x2 = -dr_x1
        dr_y2 = -dr_y1
        dr_z2 = -dr_z1

        DUMMY11 = SQRT(px1**2+py1**2+pz1**2)
        DUMMY22 = SQRT(px2**2+py2**2+pz2**2)
        DUMMY33 = DUMMY11*DUMMY22

        u1u2 = (px1*px2 + py1*py2 + pz1*pz2)/DUMMY33
        ru1 = ((x2-x1)*px1 + (y2-y1)*py1 + (z2-z1)*pz1)/(r*DUMMY11)
        ru2 = ((x2-x1)*px2 + (y2-y1)*py2 + (z2-z1)*pz2)/(r*DUMMY22)
!        rcutu1 =((x2-x1)*px1 + (y2-y1)*py1 + (z2-z1)*pz1)/(rcut*DUMMY11)  ! sigma doesn't depend on the actual distance
!        rcutu2 =((x2-x1)*px2 + (y2-y1)*py2 + (z2-z1)*pz2)/(rcut*DUMMY22) 

! get the value of sigmacut

!        DUMMY1=1.0D0/(1+khi*u1u2)
!        DUMMY2=1.0D0/(1-khi*u1u2)
!        DUMMY3=1.0D0-1.0D0/2.0D0*khi*((rcutu1+rcutu2)**2*DUMMY1+(rcutu1-rcutu2)**2*DUMMY2)
!        DUMMY4=DUMMY3
!        sigmacut= sigma0*1.0D0/(SQRT(DUMMY3))
 
!        LJ1cc=1.0D0/(rcut-sigmacut+1.0D0)
!        LJ2cc=LJ1cc**2
!        LJ3cc=LJ2cc*LJ1cc
!        LJ4cc=LJ2cc**2
!        LJ6cc=LJ4cc*LJ2cc
!        LJ7cc=LJ6cc*LJ1cc
!        LJ12cc=LJ6cc**2
!        LJ13cc=LJ12cc*LJ1cc


!((6.0D0*LJ12c - 3.0D0*LJ6c)*((r-sigma+1.0D0)/(rcut-sigmacut+1.0D0))**2 - 7.0D0*LJ12c + 4.0D0*LJ6c) 
!               (6.0D0*LJ12c - 3.0D0*LJ6c)*((r-sigma+1.0D0)/(rcut-sigmacut+1.0D0))**2 - 7.0D0*LJ12c + 4.0D0*LJ6c

!        DUMMY1=1.0D0/SQRT(1.0D0-khi**2*(u1u2)**2)
!        epsilonnu=epsilon0*DUMMY1
!        DUMMY5=DUMMY1

!        DUMMY1=1.0D0/(1.0D0 + khi1*u1u2)
!        DUMMY2=1.0D0/(1.0D0 - khi1*u1u2)
!        epsilonmu=1.0D0-khi1/2.0D0*((rcutu1+rcutu2)**2*DUMMY1 + (rcutu1-rcutu2)**2*DUMMY2)
!        epsilon=(epsilonnu**nu)*(epsilonmu**mu)   

!        LJ1=1.0D0/(rcut-sigma+1.0D0)
!        LJ2=LJ1**2
!        LJ3=LJ2*LJ1
!        LJ4=LJ2**2
!        LJ6=LJ4*LJ2
!        LJ7=LJ6*LJ1
!        LJ12=LJ6**2
!        LJ13=LJ12*LJ1

!        FLJ=LJ12-LJ6

!        Vcut = epsilon*FLJ

        DUMMY1=1.0D0/(1+khi*u1u2)
        DUMMY2=1.0D0/(1-khi*u1u2)
        DUMMY3=1.0D0-1.0D0/2.0D0*khi*((ru1+ru2)**2*DUMMY1+(ru1-ru2)**2*DUMMY2)
        DUMMY4=DUMMY3
        sigma=sigma0*1.0D0/(SQRT(DUMMY3))
 
        DUMMY1=1.0D0/SQRT(1.0D0-khi**2*(u1u2)**2)
        epsilonnu=epsilon0*DUMMY1
        DUMMY5=DUMMY1

        DUMMY1=1.0D0/(1.0D0 + khi1*u1u2)
        DUMMY2=1.0D0/(1.0D0 - khi1*u1u2)
        epsilonmu=1.0D0-khi1/2.0D0*((ru1+ru2)**2*DUMMY1 + (ru1-ru2)**2*DUMMY2)
        epsilon=(epsilonnu**nu)*(epsilonmu**mu)   

! ! LJ-type terms for r = rcut
!        LJ1c=1.0D0/(rcut-sigma+1.0D0)
!        LJ2c=LJ1c**2
!        LJ3c=LJ2c*LJ1c
!        LJ4c=LJ2c**2
!        LJ6c=LJ4c*LJ2c
!        LJ7c=LJ6c*LJ1c
!        LJ12c=LJ6c**2
!        LJ13c=LJ12c*LJ1c

! now for the original Gay-Berne potential
        LJ1=1.0D0/(r-sigma+1.0D0)
        LJ2=LJ1**2
        LJ3=LJ2*LJ1
        LJ4=LJ2**2
        LJ6=LJ4*LJ2
        LJ7=LJ6*LJ1
        LJ12=LJ6**2
        LJ13=LJ12*LJ1

        FLJ=LJ12-LJ6

!   ! Stoddard-Ford type cutoff part of the potential (sigma has to be calculated first!!! damn)
       
!        FLJc =(6.0D0*LJ12c - 3.0D0*LJ6c)*((r-sigma+1.0D0)/(rcut-sigma+1.0D0))**2 - 7.0D0*LJ12c + 4.0D0*LJ6c + sigma*LJ7cc*(2.0D0*LJ6cc-1.0D0)
!        FLJshift = sigma0*LJ7cc*(2.0D0*LJ6cc-1.0D0)

!         WRITE(*,*) 'FLJc:', r, sigma, rcut, sigmacut 

        VDUMMY=epsilon*FLJ      
!        VDUMMY = ((r-sigma+1.0D0)/(rcut-sigma+1.0D0))**2 
!        IF(VDUMMY>1.0D5) THEN
!        VDUMMY = -1.0D20
!        V(:) = 1.0D0
!        GOTO 999
!        END IF
 !       WRITE(*,*) 'sf344> energy debug',  DUMMY3

    ! derivatives follow 

        DUMMY1 = DUMMY11**3
        DUMMY2 = DUMMY22**3
        DUMMY3 = (px1*px2 + py1*py2 + pz1*pz2)

        du1u2_px1 = px2/DUMMY33 - px1*DUMMY3/(DUMMY1*DUMMY22)
        du1u2_py1 = py2/DUMMY33 - py1*DUMMY3/(DUMMY1*DUMMY22)
        du1u2_pz1 = pz2/DUMMY33 - pz1*DUMMY3/(DUMMY1*DUMMY22)
        du1u2_px2 = px1/DUMMY33 - px2*DUMMY3/(DUMMY2*DUMMY11)
        du1u2_py2 = py1/DUMMY33 - py2*DUMMY3/(DUMMY2*DUMMY11)
        du1u2_pz2 = pz1/DUMMY33 - pz2*DUMMY3/(DUMMY2*DUMMY11)

        DUMMY1 = ((x2-x1)*px1 + (y2-y1)*py1 + (z2-z1)*pz1)
        DUMMY2 = ((x2-x1)*px2 + (y2-y1)*py2 + (z2-z1)*pz2)
        r2 = r**2

        dru1_x1 = (-px1*r - dr_x1*DUMMY1)/(r2*DUMMY11)
        dru1_x2 = ( px1*r - dr_x2*DUMMY1)/(r2*DUMMY11)
        dru1_y1 = (-py1*r - dr_y1*DUMMY1)/(r2*DUMMY11)
        dru1_y2 = ( py1*r - dr_y2*DUMMY1)/(r2*DUMMY11)
        dru1_z1 = (-pz1*r - dr_z1*DUMMY1)/(r2*DUMMY11)
        dru1_z2 = ( pz1*r - dr_z2*DUMMY1)/(r2*DUMMY11)
        dru2_x1 = (-px2*r - dr_x1*DUMMY2)/(r2*DUMMY22)
        dru2_x2 = ( px2*r - dr_x2*DUMMY2)/(r2*DUMMY22)
        dru2_y1 = (-py2*r - dr_y1*DUMMY2)/(r2*DUMMY22)
        dru2_y2 = ( py2*r - dr_y2*DUMMY2)/(r2*DUMMY22)
        dru2_z1 = (-pz2*r - dr_z1*DUMMY2)/(r2*DUMMY22)
        dru2_z2 = ( pz2*r - dr_z2*DUMMY2)/(r2*DUMMY22)


        dru1_px1 = (x2-x1)/(r*DUMMY11) - px1*((x2-x1)*px1 + (y2-y1)*py1 + (z2-z1)*pz1)/(DUMMY11**3*r)
        dru1_py1 = (y2-y1)/(r*DUMMY11) - py1*((x2-x1)*px1 + (y2-y1)*py1 + (z2-z1)*pz1)/(DUMMY11**3*r)
        dru1_pz1 = (z2-z1)/(r*DUMMY11) - pz1*((x2-x1)*px1 + (y2-y1)*py1 + (z2-z1)*pz1)/(DUMMY11**3*r)
        dru2_px2 = (x2-x1)/(r*DUMMY22) - px2*((x2-x1)*px2 + (y2-y1)*py2 + (z2-z1)*pz2)/(DUMMY22**3*r)
        dru2_py2 = (y2-y1)/(r*DUMMY22) - py2*((x2-x1)*px2 + (y2-y1)*py2 + (z2-z1)*pz2)/(DUMMY22**3*r)
        dru2_pz2 = (z2-z1)/(r*DUMMY22) - pz2*((x2-x1)*px2 + (y2-y1)*py2 + (z2-z1)*pz2)/(DUMMY22**3*r)


        DUMMY1=1.0D0/(1.0D0 + khi*u1u2)
        DUMMY2=1.0D0/(1.0D0 - khi*u1u2)
        DUMMY11 = DUMMY1**2
        DUMMY22 = DUMMY2**2
        
        dsigma_px1=1.0D0/4.0D0*sigma0*khi*1.0D0/SQRT(DUMMY4**3)*((2.0D0*(ru1+ru2)*dru1_px1*(1.0D0+khi*u1u2) &
   &                 -khi*du1u2_px1*(ru1+ru2)**2)*DUMMY11 + &
          & (2.0D0*(ru1-ru2)*dru1_px1*(1.0D0-khi*u1u2)+khi*du1u2_px1*(ru1-ru2)**2)*DUMMY22)
        dsigma_py1=1.0D0/4.0D0*sigma0*khi*1.0D0/SQRT(DUMMY4**3)*((2.0D0*(ru1+ru2)*dru1_py1*(1.0D0+khi*u1u2) &
   &                 -khi*du1u2_py1*(ru1+ru2)**2)*DUMMY11 + &
          & (2.0D0*(ru1-ru2)*dru1_py1*(1.0D0-khi*u1u2)+khi*du1u2_py1*(ru1-ru2)**2)*DUMMY22)
        dsigma_pz1=1.0D0/4.0D0*sigma0*khi*1.0D0/SQRT(DUMMY4**3)*((2.0D0*(ru1+ru2)*dru1_pz1*(1.0D0+khi*u1u2) &
   &                 -khi*du1u2_pz1*(ru1+ru2)**2)*DUMMY11 + &
          & (2.0D0*(ru1-ru2)*dru1_pz1*(1.0D0-khi*u1u2)+khi*du1u2_pz1*(ru1-ru2)**2)*DUMMY22)
 


       dsigma_px2=1.0D0/4.0D0*sigma0*khi*1.0D0/SQRT(DUMMY4**3)*((2.0D0*(ru1+ru2)*dru2_px2*(1.0D0+khi*u1u2) &
   &                               -khi*du1u2_px2*(ru1+ru2)**2)*DUMMY11 + &
      & (2.0D0*(ru1-ru2)*(-dru2_px2)*(1.0D0-khi*u1u2)+khi*du1u2_px2*(ru1-ru2)**2)*DUMMY22)
       dsigma_py2=1.0D0/4.0D0*sigma0*khi*1.0D0/SQRT(DUMMY4**3)*((2.0D0*(ru1+ru2)*dru2_py2*(1.0D0+khi*u1u2) &
   &                               -khi*du1u2_py2*(ru1+ru2)**2)*DUMMY11 + &
      & (2.0D0*(ru1-ru2)*(-dru2_py2)*(1.0D0-khi*u1u2)+khi*du1u2_py2*(ru1-ru2)**2)*DUMMY22)
       dsigma_pz2=1.0D0/4.0D0*sigma0*khi*1.0D0/SQRT(DUMMY4**3)*((2.0D0*(ru1+ru2)*dru2_pz2*(1.0D0+khi*u1u2) &
   &                               -khi*du1u2_pz2*(ru1+ru2)**2)*DUMMY11 + &
      & (2.0D0*(ru1-ru2)*(-dru2_pz2)*(1.0D0-khi*u1u2)+khi*du1u2_pz2*(ru1-ru2)**2)*DUMMY22)


        dsigma_x1=1.0D0/4.0D0*sigma0*khi*1.0D0/SQRT(DUMMY4**3)*((2.0D0*(ru1+ru2)*(dru1_x1+dru2_x1))*DUMMY1 + &
             & (2.0D0*(ru1-ru2)*(dru1_x1-dru2_x1))*DUMMY2)
        dsigma_x2=1.0D0/4.0D0*sigma0*khi*1.0D0/SQRT(DUMMY4**3)*((2.0D0*(ru1+ru2)*(dru1_x2+dru2_x2))*DUMMY1 + &
             & (2.0D0*(ru1-ru2)*(dru1_x2-dru2_x2))*DUMMY2)
        dsigma_y1=1.0D0/4.0D0*sigma0*khi*1.0D0/SQRT(DUMMY4**3)*((2.0D0*(ru1+ru2)*(dru1_y1+dru2_y1))*DUMMY1 + &
             & (2.0D0*(ru1-ru2)*(dru1_y1-dru2_y1))*DUMMY2)
        dsigma_y2=1.0D0/4.0D0*sigma0*khi*1.0D0/SQRT(DUMMY4**3)*((2.0D0*(ru1+ru2)*(dru1_y2+dru2_y2))*DUMMY1 + &
             & (2.0D0*(ru1-ru2)*(dru1_y2-dru2_y2))*DUMMY2)
        dsigma_z1=1.0D0/4.0D0*sigma0*khi*1.0D0/SQRT(DUMMY4**3)*((2.0D0*(ru1+ru2)*(dru1_z1+dru2_z1))*DUMMY1 + &
             & (2.0D0*(ru1-ru2)*(dru1_z1-dru2_z1))*DUMMY2)
        dsigma_z2=1.0D0/4.0D0*sigma0*khi*1.0D0/SQRT(DUMMY4**3)*((2.0D0*(ru1+ru2)*(dru1_z2+dru2_z2))*DUMMY1 + &
             & (2.0D0*(ru1-ru2)*(dru1_z2-dru2_z2))*DUMMY2)


         DUMMY1=1.0D0/(1.0D0 + khi1*u1u2)
         DUMMY2=1.0D0/(1.0D0 - khi1*u1u2)
         DUMMY11 = DUMMY1**2
         DUMMY22 = DUMMY2**2
    
         depsilonmu_px1=-khi1/2.0D0*((2.0D0*(ru1+ru2)*dru1_px1*(1.0D0+khi1*u1u2)-khi1*du1u2_px1*(ru1+ru2)**2)*DUMMY11 + &
   &(2.0D0*(ru1-ru2)*dru1_px1*(1.0D0-khi1*u1u2)+khi1*du1u2_px1*(ru1-ru2)**2)*DUMMY22)
         depsilonmu_py1=-khi1/2.0D0*((2.0D0*(ru1+ru2)*dru1_py1*(1.0D0+khi1*u1u2)-khi1*du1u2_py1*(ru1+ru2)**2)*DUMMY11 + &
   &(2.0D0*(ru1-ru2)*dru1_py1*(1.0D0-khi1*u1u2)+khi1*du1u2_py1*(ru1-ru2)**2)*DUMMY22)
         depsilonmu_pz1=-khi1/2.0D0*((2.0D0*(ru1+ru2)*dru1_pz1*(1.0D0+khi1*u1u2)-khi1*du1u2_pz1*(ru1+ru2)**2)*DUMMY11 + &
   &(2.0D0*(ru1-ru2)*dru1_pz1*(1.0D0-khi1*u1u2)+khi1*du1u2_pz1*(ru1-ru2)**2)*DUMMY22)

         depsilonmu_px2=-khi1/2.0D0*((2.0D0*(ru1+ru2)*dru2_px2*(1.0D0+khi1*u1u2)-khi1*du1u2_px2*(ru1+ru2)**2)*DUMMY11 + &
   &(2.0D0*(ru1-ru2)*(-dru2_px2)*(1.0D0-khi1*u1u2)+khi1*du1u2_px2*(ru1-ru2)**2)*DUMMY22)
         depsilonmu_py2=-khi1/2.0D0*((2.0D0*(ru1+ru2)*dru2_py2*(1.0D0+khi1*u1u2)-khi1*du1u2_py2*(ru1+ru2)**2)*DUMMY11 + &
   &(2.0D0*(ru1-ru2)*(-dru2_py2)*(1.0D0-khi1*u1u2)+khi1*du1u2_py2*(ru1-ru2)**2)*DUMMY22)
         depsilonmu_pz2=-khi1/2.0D0*((2.0D0*(ru1+ru2)*dru2_pz2*(1.0D0+khi1*u1u2)-khi1*du1u2_pz2*(ru1+ru2)**2)*DUMMY11 + &
   &(2.0D0*(ru1-ru2)*(-dru2_pz2)*(1.0D0-khi1*u1u2)+khi1*du1u2_pz2*(ru1-ru2)**2)*DUMMY22)
 
         depsilonmu_x1=-khi1/2.0D0*(2.0D0*(ru1+ru2)*(dru1_x1+dru2_x1)*DUMMY1 + 2.0D0*(ru1-ru2)*(dru1_x1-dru2_x1)*DUMMY2)
         depsilonmu_x2=-khi1/2.0D0*(2.0D0*(ru1+ru2)*(dru1_x2+dru2_x2)*DUMMY1 + 2.0D0*(ru1-ru2)*(dru1_x2-dru2_x2)*DUMMY2)
         depsilonmu_y1=-khi1/2.0D0*(2.0D0*(ru1+ru2)*(dru1_y1+dru2_y1)*DUMMY1 + 2.0D0*(ru1-ru2)*(dru1_y1-dru2_y1)*DUMMY2)
         depsilonmu_y2=-khi1/2.0D0*(2.0D0*(ru1+ru2)*(dru1_y2+dru2_y2)*DUMMY1 + 2.0D0*(ru1-ru2)*(dru1_y2-dru2_y2)*DUMMY2)
         depsilonmu_z1=-khi1/2.0D0*(2.0D0*(ru1+ru2)*(dru1_z1+dru2_z1)*DUMMY1 + 2.0D0*(ru1-ru2)*(dru1_z1-dru2_z1)*DUMMY2)
         depsilonmu_z2=-khi1/2.0D0*(2.0D0*(ru1+ru2)*(dru1_z2+dru2_z2)*DUMMY1 + 2.0D0*(ru1-ru2)*(dru1_z2-dru2_z2)*DUMMY2)
   
         depsilonnu_px1=epsilon0*khi**2*DUMMY5**3*u1u2*du1u2_px1
         depsilonnu_py1=epsilon0*khi**2*DUMMY5**3*u1u2*du1u2_py1
         depsilonnu_pz1=epsilon0*khi**2*DUMMY5**3*u1u2*du1u2_pz1

         depsilonnu_px2=epsilon0*khi**2*DUMMY5**3*u1u2*du1u2_px2
         depsilonnu_py2=epsilon0*khi**2*DUMMY5**3*u1u2*du1u2_py2
         depsilonnu_pz2=epsilon0*khi**2*DUMMY5**3*u1u2*du1u2_pz2

         depsilon_px1=nu*epsilonnu**(nu-1.0D0)*depsilonnu_px1*epsilonmu**mu &
   &                       + mu*epsilonmu**(mu-1.0D0)*depsilonmu_px1*epsilonnu**nu
         depsilon_py1=nu*epsilonnu**(nu-1.0D0)*depsilonnu_py1*epsilonmu**mu & 
   &                       + mu*epsilonmu**(mu-1.0D0)*depsilonmu_py1*epsilonnu**nu
         depsilon_pz1=nu*epsilonnu**(nu-1.0D0)*depsilonnu_pz1*epsilonmu**mu & 
   &                       + mu*epsilonmu**(mu-1.0D0)*depsilonmu_pz1*epsilonnu**nu

         depsilon_px2=nu*epsilonnu**(nu-1.0D0)*depsilonnu_px2*epsilonmu**mu &
   &                       + mu*epsilonmu**(mu-1.0D0)*depsilonmu_px2*epsilonnu**nu
         depsilon_py2=nu*epsilonnu**(nu-1.0D0)*depsilonnu_py2*epsilonmu**mu &
   &                       + mu*epsilonmu**(mu-1.0D0)*depsilonmu_py2*epsilonnu**nu
         depsilon_pz2=nu*epsilonnu**(nu-1.0D0)*depsilonnu_pz2*epsilonmu**mu &
   &                       + mu*epsilonmu**(mu-1.0D0)*depsilonmu_pz2*epsilonnu**nu


         depsilon_x1= mu*epsilonmu**(mu-1.0D0)*depsilonmu_x1*epsilonnu**nu
         depsilon_x2= mu*epsilonmu**(mu-1.0D0)*depsilonmu_x2*epsilonnu**nu
         depsilon_y1= mu*epsilonmu**(mu-1.0D0)*depsilonmu_y1*epsilonnu**nu
         depsilon_y2= mu*epsilonmu**(mu-1.0D0)*depsilonmu_y2*epsilonnu**nu
         depsilon_z1= mu*epsilonmu**(mu-1.0D0)*depsilonmu_z1*epsilonnu**nu
         depsilon_z2= mu*epsilonmu**(mu-1.0D0)*depsilonmu_z2*epsilonnu**nu

!        ! for the Stoddard-Ford cutoff part

!        dFLJc_x1 =   2.0D0*(r-sigma+1.0D0)*((-1.0D0/r*(x2-x1)-dsigma_x1)*(rcut-sigma+1.0D0)+dsigma_x1*(r-sigma+1.0D0))*&
!                        & (1.0D0/(rcut-sigma+1.0D0))**3*(6.0D0*LJ12c - 3.0D0*LJ6c)-6.0D0*dsigma_x1*LJ7c*(2.0D0*LJ6c-1.0D0)
!        dFLJc_x2 =   2.0D0*(r-sigma+1.0D0)*(( 1.0D0/r*(x2-x1)-dsigma_x2)*(rcut-sigma+1.0D0)+dsigma_x2*(r-sigma+1.0D0))*&
!                        & (1.0D0/(rcut-sigma+1.0D0))**3*(6.0D0*LJ12c - 3.0D0*LJ6c)-6.0D0*dsigma_x2*LJ7c*(2.0D0*LJ6c-1.0D0)
!        dFLJc_y1 =   2.0D0*(r-sigma+1.0D0)*((-1.0D0/r*(y2-y1)-dsigma_y1)*(rcut-sigma+1.0D0)+dsigma_y1*(r-sigma+1.0D0))*&
!                        & (1.0D0/(rcut-sigma+1.0D0))**3*(6.0D0*LJ12c - 3.0D0*LJ6c)-6.0D0*dsigma_y1*LJ7c*(2.0D0*LJ6c-1.0D0)
!        dFLJc_y2 =   2.0D0*(r-sigma+1.0D0)*(( 1.0D0/r*(y2-y1)-dsigma_y2)*(rcut-sigma+1.0D0)+dsigma_y2*(r-sigma+1.0D0))*&
!                        & (1.0D0/(rcut-sigma+1.0D0))**3*(6.0D0*LJ12c - 3.0D0*LJ6c)-6.0D0*dsigma_y2*LJ7c*(2.0D0*LJ6c-1.0D0)
!        dFLJc_z1 =   2.0D0*(r-sigma+1.0D0)*((-1.0D0/r*(z2-z1)-dsigma_z1)*(rcut-sigma+1.0D0)+dsigma_z1*(r-sigma+1.0D0))*&
!                        & (1.0D0/(rcut-sigma+1.0D0))**3*(6.0D0*LJ12c - 3.0D0*LJ6c)-6.0D0*dsigma_z1*LJ7c*(2.0D0*LJ6c-1.0D0)
!        dFLJc_z2 =   2.0D0*(r-sigma+1.0D0)*(( 1.0D0/r*(z2-z1)-dsigma_z2)*(rcut-sigma+1.0D0)+dsigma_z2*(r-sigma+1.0D0))*&
!                        & (1.0D0/(rcut-sigma+1.0D0))**3*(6.0D0*LJ12c - 3.0D0*LJ6c)-6.0D0*dsigma_z2*LJ7c*(2.0D0*LJ6c-1.0D0)
 
!        dFLJc_px1 = -2.0D0*(r-sigma+1.0D0)*(-dsigma_px1*(rcut-sigma+1.0D0)+dsigma_px1*(r-sigma+1))*&
!                        & (1.0D0/(rcut-sigma+1.0D0))**3*(6.0D0*LJ12c-3.0D0*LJ6c)-6.0D0*dsigma_px1*LJ7c*(2.0D0*LJ6c-1.0D0) 
!        dFLJc_py1 = -2.0D0*(r-sigma+1.0D0)*(-dsigma_py1*(rcut-sigma+1.0D0)+dsigma_py1*(r-sigma+1))*&
!                        & (1.0D0/(rcut-sigma+1.0D0))**3*(6.0D0*LJ12c-3.0D0*LJ6c)-6.0D0*dsigma_py1*LJ7c*(2.0D0*LJ6c-1.0D0)
!        dFLJc_pz1 = -2.0D0*(r-sigma+1.0D0)*(-dsigma_pz1*(rcut-sigma+1.0D0)+dsigma_pz1*(r-sigma+1))*&
!                        & (1.0D0/(rcut-sigma+1.0D0))**3*(6.0D0*LJ12c-3.0D0*LJ6c)-6.0D0*dsigma_pz1*LJ7c*(2.0D0*LJ6c-1.0D0)
!        dFLJc_px2 = -2.0D0*(r-sigma+1.0D0)*(-dsigma_px2*(rcut-sigma+1.0D0)+dsigma_px2*(r-sigma+1))*&
!                        & (1.0D0/(rcut-sigma+1.0D0))**3*(6.0D0*LJ12c-3.0D0*LJ6c)-6.0D0*dsigma_px2*LJ7c*(2.0D0*LJ6c-1.0D0)
!        dFLJc_py2 = -2.0D0*(r-sigma+1.0D0)*(-dsigma_py2*(rcut-sigma+1.0D0)+dsigma_py2*(r-sigma+1))*&
!                        & (1.0D0/(rcut-sigma+1.0D0))**3*(6.0D0*LJ12c-3.0D0*LJ6c)-6.0D0*dsigma_py2*LJ7c*(2.0D0*LJ6c-1.0D0)
!        dFLJc_pz2 = -2.0D0*(r-sigma+1.0D0)*(-dsigma_pz2*(rcut-sigma+1.0D0)+dsigma_pz2*(r-sigma+1))*&
!                        & (1.0D0/(rcut-sigma+1.0D0))**3*(6.0D0*LJ12c-3.0D0*LJ6c)-6.0D0*dsigma_pz2*LJ7c*(2.0D0*LJ6c-1.0D0)

        ! now for the original GB potential + the cutoff part

         dFLJ_px1=6.0D0*dsigma_px1*LJ7*(2.0D0*LJ6-1.0D0)! + dFLJc_px1
         dFLJ_py1=6.0D0*dsigma_py1*LJ7*(2.0D0*LJ6-1.0D0)! + dFLJc_py1
         dFLJ_pz1=6.0D0*dsigma_pz1*LJ7*(2.0D0*LJ6-1.0D0)! + dFLJc_pz1
         dFLJ_px2=6.0D0*dsigma_px2*LJ7*(2.0D0*LJ6-1.0D0)! + dFLJc_px2
         dFLJ_py2=6.0D0*dsigma_py2*LJ7*(2.0D0*LJ6-1.0D0)! + dFLJc_py2
         dFLJ_pz2=6.0D0*dsigma_pz2*LJ7*(2.0D0*LJ6-1.0D0)! + dFLJc_pz2

         dFLJ_x1=-6.0D0*LJ7*(2.0D0*LJ6-1.0D0)*(-1.0D0/r*(x2-x1)-dsigma_x1)! + dFLJc_x1
         dFLJ_x2=-6.0D0*LJ7*(2.0D0*LJ6-1.0D0)*( 1.0D0/r*(x2-x1)-dsigma_x2)! + dFLJc_x2
         dFLJ_y1=-6.0D0*LJ7*(2.0D0*LJ6-1.0D0)*(-1.0D0/r*(y2-y1)-dsigma_y1)! + dFLJc_y1
         dFLJ_y2=-6.0D0*LJ7*(2.0D0*LJ6-1.0D0)*( 1.0D0/r*(y2-y1)-dsigma_y2)! + dFLJc_y2
         dFLJ_z1=-6.0D0*LJ7*(2.0D0*LJ6-1.0D0)*(-1.0D0/r*(z2-z1)-dsigma_z1)! + dFLJc_z1
         dFLJ_z2=-6.0D0*LJ7*(2.0D0*LJ6-1.0D0)*( 1.0D0/r*(z2-z1)-dsigma_z2)! + dFLJc_z2


         dV_x1=depsilon_x1*(FLJ) + epsilon*dFLJ_x1 
         dV_x2=depsilon_x2*(FLJ) + epsilon*dFLJ_x2 
         dV_y1=depsilon_y1*(FLJ) + epsilon*dFLJ_y1 
         dV_y2=depsilon_y2*(FLJ) + epsilon*dFLJ_y2 
         dV_z1=depsilon_z1*(FLJ) + epsilon*dFLJ_z1 
         dV_z2=depsilon_z2*(FLJ) + epsilon*dFLJ_z2 

         dV_px1=depsilon_px1*(FLJ) + epsilon*dFLJ_px1 
         dV_py1=depsilon_py1*(FLJ) + epsilon*dFLJ_py1 
         dV_pz1=depsilon_pz1*(FLJ) + epsilon*dFLJ_pz1 

         dV_px2=depsilon_px2*(FLJ) + epsilon*dFLJ_px2 
         dV_py2=depsilon_py2*(FLJ) + epsilon*dFLJ_py2 
         dV_pz2=depsilon_pz2*(FLJ) + epsilon*dFLJ_pz2 

!END IF
!        WRITE(*,*) 'sf344: function values:'
!        WRITE(*,*) VDUMMY,FLJc,r
!        WRITE(*,*) LJ12-LJ6+((6.0D0*LJ12c - 3.0D0*LJ6c)*((r-sigma+1.0D0)/(rcut-sigmacut+1.0D0))**2)- 7.0D0*LJ12c + 4.0D0*LJ6c
!        WRITE(*,*) FLJ +((6.0D0*LJ12c - 3.0D0*LJ6c)*((r-sigma+1.0D0)/(rcut-sigmacut+1.0D0))**2- 7.0D0*LJ12c + 4.0D0*LJ6c)
!        WRITE(*,*) sigmacut, sigma
!        WRITE(*,*) dV_px1, dV_px2, dV_py1, dV_py2, VDUMMY 
!        STOP

   
         EGB=EGB+VDUMMY
    

!        VT(j) = VT(j) + VDUMMY
!        VT(k) = VT(k) + VDUMMY        ! pair potentials

         V(J2-2)=V(J2-2)+dV_x1
         V(K2-2)=V(K2-2)+dV_x2
         V(J2-1)=V(J2-1)+dV_y1
         V(K2-1)=V(K2-1)+dV_y2
         V(J2  )=V(J2  )+dV_z1
         V(K2  )=V(K2  )+dV_z2
         V(3*REALNATOMS+J2-2)=V(3*REALNATOMS+J2-2)+dV_px1
         V(3*REALNATOMS+K2-2)=V(3*REALNATOMS+K2-2)+dV_px2
         V(3*REALNATOMS+J2-1)=V(3*REALNATOMS+J2-1)+dV_py1
         V(3*REALNATOMS+K2-1)=V(3*REALNATOMS+K2-1)+dV_py2
         V(3*REALNATOMS+J2  )=V(3*REALNATOMS+J2  )+dV_pz1
         V(3*REALNATOMS+K2  )=V(3*REALNATOMS+K2  )+dV_pz2
!WRITE(*,*) 'sf344> ', V(:) 

!        NUMDER(i+1) = ru2
       
! END DO
!        DERIV = (NUMDER(3)-NUMDER(1))/0.0002
!        WRITE(*,'(A, 3F18.10)') 'sf344> numerical derivatives', DERIV, dru2_py2, DERIV-dru2_py2
!        WRITE(*,*) X(:)
!        STOP
    END DO
END DO
!WRITE(*,*) 'OVERLAPTEST:', OVERLAPTEST
!WRITE(*,*) 'sf344> checking values', r, sigma, r-sigma
! check for overlap
!        WRITE(*,*) x1,y1,z1,px1,py1,pz1,x2,y2,z2,px2,py2,pz2,GBANISOTROPYR/2.0D0,0.5D0,0.5D0

!        ECFvalue = 1.0D0

IF(OVERLAPTEST) THEN
!WRITE(*,*) 'in here 1'
DO 
   OVERLAP = .FALSE.

   DO J1 = 1,REALNATOMS
      J2 = 3*J1

      DO K1 = J1+1,REALNATOMS
         K2 = 3*K1
         r = SQRT((X(J2-2)-X(K2-2))**2+(X(J2-1)-X(K2-1))**2+(X(J2)-X(K2))**2)
         IF (r<GBANISOTROPYR) THEN
          DO
           CALL ECF(OVERLAPTEST,ECFvalue,X(J2-2),X(K2-2),X(J2-1),X(K2-1),X(J2),X(K2),&
                        & X(3*REALNATOMS+J2-2),X(3*REALNATOMS+K2-2),X(3*REALNATOMS+J2-1),X(3*REALNATOMS+K2-1),&
                        & X(3*REALNATOMS+J2),X(3*REALNATOMS+K2),GBANISOTROPYR/2.0D0,0.5D0,0.5D0)
!            WRITE(*,*) 'in here 2', OVERLAPTEST, ECFvalue
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


SUBROUTINE OLDECF(OVERLAPTEST,ECFvalue,x1,x2,y1,y2,z1,z2,alpha1,alpha2,beta1,beta2,a,b)


!Calculates the value of the ellipsoid contact function
use commons
!use keyword 
IMPLICIT NONE

DOUBLE PRECISION       ::x1,x2,y1,y2,z1,z2,alpha1,alpha2,beta1,beta2,a,b
DOUBLE PRECISION  ::lambda,C_alpha1,C_beta1,C_alpha2,C_beta2,S_alpha1,S_beta1,S_alpha2,S_beta2,A1(3,3),B1(3,3),C1(3,3),D1(6),&
&    Rvect(3),RvectT(1,3),S_lambda(1),ECFvalue
DOUBLE PRECISION  ::S_lambdamax(101,1)
INTEGER    ::i,j,K
LOGICAL    ::OVERLAPTEST

S_lambdamax(:,1)=0.0D0

a=GBANISOTROPYR/2.0D0
b=0.5D0


S_alpha1=SIN(alpha1)
S_beta1=SIN(beta1)
S_alpha2=SIN(alpha2)
S_beta2=SIN(beta2)
C_alpha1=COS(alpha1)
C_beta1=COS(beta1)
C_alpha2=COS(alpha2)
C_beta2=COS(beta2)


A1(1,1) = ((C_alpha1*C_beta1)**2)/a**2 +((S_beta1)**2)/b**2 +((C_beta1*S_alpha1)**2)/b**2
A1(1,2) = 0
A1(1,3) = 0
A1(2,1) = A1(1,2)
A1(2,2) = ((C_alpha1*S_beta1)**2)/a**2 +((C_beta1)**2)/b**2 +((S_alpha1*S_beta1)**2)/b**2
A1(2,3) = 0
A1(3,1) = A1(1,3)
A1(3,2) = A1(2,3)
A1(3,3) = (S_alpha1)**2/a**2 +(C_alpha1)**2/b**2

B1(1,1) = ((C_alpha2*C_beta2)**2)/a**2 +((S_beta2)**2)/b**2 +((C_beta2*S_alpha2)**2)/b**2
B1(1,2) = 0
B1(1,3) = 0
B1(2,1) = B1(1,2)
B1(2,2) = ((C_alpha2*S_beta2)**2)/a**2 +((C_beta2)**2)/b**2 +((S_alpha2*S_beta2)**2)/b**2
B1(2,3) = 0
B1(3,1) = B1(1,3)
B1(3,2) = B1(2,3)
B1(3,3) = (S_alpha2)**2/a**2 +(C_alpha2)**2/b**2

!diagonal matrix inversion



DO i=1,3
    A1(i,i)=1/A1(i,i)
    B1(i,i)=1/B1(i,i)
END DO



!WRITE(*,*) A1(:,:)
!WRITE(*,*)
!WRITE(*,*) B1(:,:)


Rvect(1)=x2-x1
Rvect(2)=y2-y1
Rvect(3)=z2-z1
RvectT(1,1)=Rvect(1)
RvectT(1,2)=Rvect(2)
RvectT(1,3)=Rvect(3)

DO k=0,100
    lambda=0.01*k
    
    DO i=1,3
 C1(i,i)=(1-lambda)*A1(i,i)+lambda*B1(i,i)
 C1(i,i)=1/C1(i,i)
    END DO
    
    S_lambda=lambda*(1-lambda)*MATMUL(RvectT,MATMUL(C1,Rvect))
    S_lambdamax(k+1,1)=S_lambda(1)
    IF(k>0.AND.S_lambdamax(k+1,1)<=S_lambdamax(k,1)) EXIT
END DO



ECFvalue=S_lambdamax(k,1)

IF(ECFvalue<1.1D0) THEN
 OVERLAPTEST=.TRUE.
ELSE
 OVERLAPTEST=.FALSE.
END IF

!WRITE(*,*) k
!WRITE(*,*) OVERLAPTEST

END SUBROUTINE OLDECF


SUBROUTINE OLDECFcheck(X,GBOVERLAP)
use commons
!use keyword
IMPLICIT NONE

DOUBLE PRECISION  :: X(3*NATOMS),x1,x2,y1,y2,z1,z2,px1,px2,py1,py2,pz1,pz2,ECFvalue,a,b,c,R
INTEGER    :: J1,J2,J3,J4,REALNATOMS
LOGICAL    :: GBOVERLAP,OVERLAPTEST

IF(GAYBERNET) THEN
        a=(GBANISOTROPYR/2.0D0)**2
        b=(0.5D0)**2
        c=(0.5D0)**2
ELSE IF(PARAMONOVT) THEN
        a=PARAMa1**2
        b=PARAMb1**2
        c=PARAMc1**2
ELSE
 WRITE(*,*) 'ECFcheck> this routine is intended to check overlap between ellipsoids, please specify &
&                which system you want to use (GAYBERNET or PARAMONOVT)'
  STOP
            
END IF

REALNATOMS=NATOMS/2

GBOVERLAP=.FALSE.


!WRITE(*,*) a,b

DO J1=1,REALNATOMS
    J2=3*J1
      
          x1=X(J2-2)
          y1=X(J2-1)
          z1=X(J2)
          px1=X(3*REALNATOMS+J2-2)
          py1=X(3*REALNATOMS+J2-1)
          pz1=X(3*REALNATOMS+J2  )

  DO J3=J1+1,REALNATOMS
      J4=3*J3
  OVERLAPTEST=.FALSE.
       x2=X(J4-2)
       y2=X(J4-1)
       z2=X(J4)
       px2=X(3*REALNATOMS+J4-2)
       py2=X(3*REALNATOMS+J4-1)
       pz2=X(3*REALNATOMS+J4  )

       r=SQRT((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
     
      CALL ECF(OVERLAPTEST,ECFvalue,x1,x2,y1,y2,z1,z2,px1,px2,py1,py2,pz1,pz2,a,b,c)
!      WRITE(*,*) ECFvalue, a,b
!      WRITE(*,*) x1,y1,z1,x2,y2,z2,alpha1,beta1
!      WRITE(*,*) alpha2,beta2
      IF (OVERLAPTEST) THEN
         GBOVERLAP=.TRUE.
!   WRITE(*,*) J1, J3
      END IF
!      WRITE(*,*) OVERLAPTEST
      END DO    
!      WRITE(*,*) GBOVERLAP
!    WRITE(*,*) X(:)
      
END DO


END SUBROUTINE OLDECFcheck

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

! numeric second derivatives for the Paramonov potential


SUBROUTINE PARAMONOVSECDER(OLDX,STEST)
use commons
!use keyword
!use modhess
implicit none

DOUBLE PRECISION  :: V(3*NATOMS),EGB,X(3*NATOMS),OLDX(3*NATOMS),VTEMP(2,3*NATOMS),ksi
INTEGER    :: i,j
LOGICAL    :: GTEST,STEST


ksi=0.0001D0
X(:)=OLDX(:)

!WRITE(*,*) 'sf344> in Paramonovsecder'
DO i=1,3*NATOMS

!       IF ((i.GT.3*NATOMS/2).and.(MOD(i,3).EQ.0)) THEN
!        VTEMP(:,:) = 0.0D0
!       ELSE
         X(i)=X(i)-ksi
 
         CALL OLDPARAMONOV(X,V,EGB,GTEST,.FALSE.)
 
         VTEMP(1,:)=V(:)
 
         X(i)=X(i)+2.0D0*ksi

         CALL OLDPARAMONOV(X,V,EGB,GTEST,.FALSE.)
 
         VTEMP(2,:)=V(:)

                DO j=1,3*NATOMS
                        HESS(j,i)=(VTEMP(2,j)-VTEMP(1,j))/(2.0D0*ksi)
                END DO
!        END IF
END DO
!WRITE(*,*) 'sf344> exiting Paramonovsecder'
!WRITE(*,*) 'HESSIAN:'
!WRITE(*,'12F10.3') HESS(:,:)

END SUBROUTINE PARAMONOVSECDER

SUBROUTINE GAYBERNESECDER(OLDX,STEST)
use commons
!use keywords
!use modhess
implicit none

DOUBLE PRECISION  :: V(3*NATOMS),EGB,X(3*NATOMS),OLDX(3*NATOMS),VTEMP(2,3*NATOMS),ksi
INTEGER    :: i,j
LOGICAL    :: GTEST,STEST


ksi=0.000001D0
X(:)=OLDX(:)

!WRITE(*,*) 'sf344> in GayBerneSecDer'
DO i=1,3*NATOMS

         X(i)=X(i)-ksi
 
         CALL GAYBERNE(X,V,EGB,GTEST,.FALSE.)
 
         VTEMP(1,:)=V(:)
 
         X(i)=X(i)+2.0D0*ksi

         CALL GAYBERNE(X,V,EGB,GTEST,.FALSE.)
 
         VTEMP(2,:)=V(:)

                DO j=1,3*NATOMS
                        HESS(j,i)=(VTEMP(2,j)-VTEMP(1,j))/(2.0D0*ksi)
                END DO
END DO
!WRITE(*,*) 'sf344> exiting GayBerneSecDer'
!WRITE(*,*) 'HESSIAN:'
!WRITE(*,'12F10.3') HESS(:,:)

END SUBROUTINE GAYBERNESECDER


SUBROUTINE TAKESTEPELLIPSOIDS(NP)
! take step routine for anisotropic particles


USE commons
IMPLICIT NONE

INTEGER                 :: NP, J1, J2, JMAX, J3,K1,K2,K3, JMAX2,j,k,i, REALNATOMS,JP,CLOSESTATOMINDEX
LOGICAL                 :: OVERLAPT, OVERLAPT2, GTEST, DISSOCIATED(NATOMS/2),DISSOC
DOUBLE PRECISION        ::x1,y1,z1,x2,y2,z2,px1,px2,py1,py2,pz1,pz2,DIST(3*NATOMS/2),r,XMASS,YMASS,ZMASS,DMAX,VMAX,VMAX2
DOUBLE PRECISION        :: VMIN,CMMAX,CMDIST(NATOMS/2),CDIST(NATOMS/2),LOCALSTEP,theta,phi,psi,PI,ECFvalue,a,b,c,MINDISTANCE

DOUBLE PRECISION        :: EGB, DPRAND,RANDOM,DUMMY, DUMMY1, DUMMY2, DUMMY3,DUMMY4, DUMMY5, RANDOMZ,RVEC(NATOMS/2*(NATOMS/2-1)/2)





!a=PARAMa1**2
!b=PARAMb1**2
!c=PARAMc1**2
IF(GAYBERNET) THEN
        a=(GBANISOTROPYR/2)**2
        b=0.5D0**2
        c=0.5D0**2
ELSE IF(PARAMONOVT) THEN
        a=PARAMa2**2
        b=PARAMb2**2
        c=PARAMc2**2
END IF

PI=ATAN(1.0D0)*4
OVERLAPT=.TRUE.

REALNATOMS=NATOMS/2



    DO J1=1,REALNATOMS
        J2=3*J1
        DUMMY2=COORDS(J2-2,NP)**2+COORDS(J2-1,NP)**2+COORDS(J2,NP)**2
        IF (DUMMY2.GT.RADIUS) THEN
            WRITE(*,'(A,I5,5F20.10)') 'J1,RAD,R**2,x,y,z:', J1, RADIUS, DUMMY2, COORDS(J2-2,NP), COORDS(J2-1,NP), COORDS(J2,NP)
            PRINT*, 'initial coordinate outside container -- increase container radius'
            STOP
        END IF
    END DO


    DO J1=1,3*NATOMS
        COORDSO(J1,NP)=COORDS(J1,NP)
    END DO

    DO J1=1,NATOMS/2
        VATO(J1,NP)=VAT(J1,NP)
    END DO





!DO j=1,NATOMS
!    WRITE(*,'(A,I2,3F18.10)') "dumping coordinates", j,COORDS(3*j-2,NP), COORDS(3*j-1,NP), COORDS(3*j,NP)
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
!  Find the most weakly bound atom, JMAX, the second most weakly bound atom, JMAX2,
!  and the pair energy of the most tightly bound atom, VMIN. An angular step is
!  taken for JMAX if its pair energy is > ASTEP*VMIN putting the atom at a radius of
!  DMAX (or CMMAX from CM of the cluster).
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


i=0
OVERLAPT=.FALSE.
DO
i=i+1
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



                        IF(VAT(J1,NP).GT.(ASTEP(NP)*VMIN).AND.(J1.EQ.JMAX).AND.(.NOT.NORESET).AND.(J1.LT.NATOMS/2).AND.i==1) THEN
!                                    WRITE(*,*) "Angular move"
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
                        ! end of angular move block
                        ELSE IF(i==1) THEN
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
                        x1=COORDS(J2-2,NP)
                        y1=COORDS(J2-1,NP)
                        z1=COORDS(J2,NP)
                    px1=COORDS(3*REALNATOMS+J2-2,NP)
                    py1=COORDS(3*REALNATOMS+J2-1,NP)
                    pz1=COORDS(3*REALNATOMS+J2  ,NP)

                END IF
        
                DO k=J1+1,NATOMS
                    K2=3*k
                
                IF(k<=REALNATOMS) THEN   
                        x2=COORDS(K2-2,NP)
                        y2=COORDS(K2-1,NP)
                        z2=COORDS(K2,NP)
                        px2=COORDS(3*REALNATOMS+K2-2,NP)
                        py2=COORDS(3*REALNATOMS+K2-1,NP)
                        pz2=COORDS(3*REALNATOMS+K2  ,NP)
                       r=SQRT((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
                END IF
    
    
                IF(k<=REALNATOMS.AND.J1<=REALNATOMS.AND.r**2<MAX(a,b,c)) THEN
                  CALL ECF(.FALSE.,ECFvalue,x1,x2,y1,y2,z1,z2,px1,px2,py1,py2,pz1,pz2,a,b,c)
!                        WRITE(*,*) ECFvalue                            
                            IF(ECFvalue<1.0D0)THEN
                              OVERLAPT2=.TRUE.
                    !         WRITE(*,*) "overlapping ellipsoids",ECFvalue,r
                    !         WRITE(*,'(A,3F18.10)') "Old x,r,sigma:",COORDS(K2-2,NP),r,sigma
                    !         WRITE(*,'(A,I2,I2)') "J1,k:",J1,k
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


! now determine whether any particle has dissociated
  DO
  DISSOC=.FALSE.
   DO J1=1,natoms/2-1
    J2=3*J1
    x1=COORDS(J2-2,NP)
    y1=COORDS(J2-1,NP)
    z1=COORDS(J2,NP)
   
    DISSOCIATED(J1)=.TRUE.
    DO K1=J1+1,natoms/2
      K2=3*K1
      x2=COORDS(K2-2,NP)
      y2=COORDS(K2-1,NP)
      z2=COORDS(K2,NP)
      RVEC(K1)=SQRT((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
      IF(RVEC(K1)<PCUTOFF) DISSOCIATED(J1)=.FALSE.
!      WRITE(*,*) 'in here, cutoff=',pcutoff,J1
    END DO
    IF(DISSOCIATED(J1)) THEN
        ! determine the closest neighbor
        MINDISTANCE = HUGE(1.0D0)
        DO K1=J1+1,natoms/2
         IF(RVEC(K1)<MINDISTANCE) THEN
           MINDISTANCE=RVEC(K1)
           CLOSESTATOMINDEX=K1
         END IF
        END DO
        WRITE(*,*) 'Atom ', J1, 'is dissociated, distance to the closest atom: ', &
        & MINDISTANCE, CLOSESTATOMINDEX
        DISSOC=.TRUE.
        ! move the dissociated atom closer to the closest atom
        COORDS(J2-2,NP)=0.5D0*(COORDS(3*CLOSESTATOMINDEX-2,NP)+COORDS(J2-2,NP))
        COORDS(J2-1,NP)=0.5D0*(COORDS(3*CLOSESTATOMINDEX-1,NP)+COORDS(J2-1,NP))
        COORDS(J2  ,NP)=0.5D0*(COORDS(3*CLOSESTATOMINDEX  ,NP)+COORDS(J2  ,NP))
    END IF
   END DO
  IF(.NOT.DISSOC) EXIT
  END DO
 IF (.NOT.OVERLAPT.AND..NOT.DISSOC) EXIT
END DO                    
!WRITE(*,*) i,"takestep moves"


END SUBROUTINE TAKESTEPELLIPSOIDS

SUBROUTINE EllipsoidsAAtoPolar(px1,py1,pz1,alpha,beta,gamma,alphadeg,betadeg,gammadeg)
! converts angle-axis coordinates px, py, pz to "polar-like" angles alpha and beta
! px=cos(alpha)*cos(alpha)
! py=cos(alpha)*sin(beta)
! pz=sin(alpha)
!use commons
implicit none

DOUBLE PRECISION        :: px1, py1, pz1,px,py,pz

DOUBLE PRECISION,intent(out)       :: alpha, beta, gamma, alphadeg, betadeg, gammadeg
DOUBLE PRECISION                   :: PI,twopi

! alpha: angle of the vector with the xy plane
! beta: angle of the vector with the x axis
! gamma: angle of rotation about the z axis (angle-axis convention),
!        provided by the magnitude of the vector (sqrt(px**2+py**2+pz**2)-1)
!        so that a vector of magnitude 1 means no rotation is performed

!IF(pz*py*pz<0) THEN
!       px = -px
!       py = -py
!       pz = -pz
!END IF

PI = atan(1.0)*4.0D0
twopi = 2.0D0*PI

     gamma = sqrt(px1**2+py1**2+pz1**2)

     px = px1/gamma
     py = py1/gamma
     pz = pz1/gamma

     IF(px.eq.0.0d0) THEN
       IF(py>=0.0D0) THEN
        alpha = PI/2           ! Euler angle alpha 
       ELSE
        alpha = -1.0D0*PI/2
       END IF
     ELSE
        IF(py>=0.0D0) THEN
          IF(px>0.0D0) THEN
            alpha = 1.0D0*atan(py/px)   ! first quadrant
          ELSE    ! px<0
            alpha = 1.0D0*atan(py/px)+PI       ! should be in the second quadrant
          END IF
        ELSE IF(py<0.0D0) THEN
          IF(px>0.0D0) THEN
            alpha = 1.0D0*atan(py/px)    ! fourth quadrant
          ELSE    ! px<0
            alpha = 1.0D0*atan(py/px)-PI            ! third quadrant
          END IF
        END IF 
     END IF
        beta = 1.0D0*acos(pz)

alphadeg=alpha*180/Pi
betadeg=beta*180/Pi
gammadeg=gamma*180/Pi

!write(*,*) 'exiting EllipsoidsAAtoPolar'
END SUBROUTINE EllipsoidsAAtoPolar


SUBROUTINE TAKESTEPPARAM(NP)
! take step routine for anisotropic particles


USE commons
IMPLICIT NONE

INTEGER                        :: NP, J1, J2, JMAX, J3,K2,K3, JMAX2,j,k,i, REALNATOMS,JP
LOGICAL                        :: OVERLAPT, OVERLAPT2, GTEST

DOUBLE PRECISION        ::x1, y1, z1, x2, y2, z2, px1,py1,pz1,px2,py2,pz2,DIST(3*NATOMS/2),r,XMASS,YMASS,ZMASS,DMAX,VMAX,VMAX2
DOUBLE PRECISION        :: VMIN,CMMAX,CMDIST(NATOMS/2),CDIST(NATOMS/2),LOCALSTEP,theta,phi,psi,PI,ECFvalue,a,b,c

DOUBLE PRECISION        :: EGB, DPRAND,RANDOM,DUMMY, DUMMY1, DUMMY2, DUMMY3,DUMMY4, DUMMY5, RANDOMZ




!a=PARAMa1**2
!b=PARAMb1**2
!c=PARAMc1**2

a=PARAMa2**2
b=PARAMb2**2
c=PARAMc2**2


PI=ATAN(1.0D0)*4
OVERLAPT=.TRUE.

REALNATOMS=NATOMS/2



    DO J1=1,REALNATOMS
        J2=3*J1
        DUMMY2=COORDS(J2-2,NP)**2+COORDS(J2-1,NP)**2+COORDS(J2,NP)**2
        IF (DUMMY2.GT.RADIUS) THEN
            WRITE(*,'(A,I5,5F20.10)') 'J1,RAD,R**2,x,y,z:', J1, RADIUS, DUMMY2, COORDS(J2-2,NP), COORDS(J2-1,NP), COORDS(J2,NP)
            PRINT*, 'initial coordinate outside container -- increase container radius'
            STOP
        END IF
    END DO


    DO J1=1,3*NATOMS
        COORDSO(J1,NP)=COORDS(J1,NP)
    END DO

    DO J1=1,NATOMS/2
        VATO(J1,NP)=VAT(J1,NP)
    END DO





!DO j=1,NATOMS
!    WRITE(*,'(A,I2,3F18.10)') "dumping coordinates", j,COORDS(3*j-2,NP), COORDS(3*j-1,NP), COORDS(3*j,NP)
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
!  Find the most weakly bound atom, JMAX, the second most weakly bound atom, JMAX2,
!  and the pair energy of the most tightly bound atom, VMIN. An angular step is
!  taken for JMAX if its pair energy is > ASTEP*VMIN putting the atom at a radius of
!  DMAX (or CMMAX from CM of the cluster).
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




i=0
OVERLAPT=.FALSE.
DO
i=i+1
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



                        IF(VAT(J1,NP).GT.(ASTEP(NP)*VMIN).AND.(J1.EQ.JMAX).AND.(.NOT.NORESET).AND.(J1.LT.NATOMS/2).AND.i==1) THEN
!                                    WRITE(*,*) "Angular move"
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
                        ! end of angular move block
                        ELSE IF(i==1) THEN
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
                        x1=COORDS(J2-2,NP)
                        y1=COORDS(J2-1,NP)
                        z1=COORDS(J2,NP)
                    px1=COORDS(3*REALNATOMS+J2-2,NP)
                    py1=COORDS(3*REALNATOMS+J2-1,NP)
                    pz1=COORDS(3*REALNATOMS+J2  ,NP)
                END IF
        
                DO k=J1+1,NATOMS
                    K2=3*k
                
                    IF(k<=REALNATOMS) THEN    
                            x2=COORDS(K2-2,NP)
                        y2=COORDS(K2-1,NP)
                        z2=COORDS(K2,NP)
                        px2=COORDS(3*REALNATOMS+K2-2,NP)
                        py2=COORDS(3*REALNATOMS+K2-1,NP)
                        pz2=COORDS(3*REALNATOMS+K2  ,NP)

                            r=DSQRT((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
                    END IF
                    
                    
                    IF(k<=REALNATOMS.AND.J1<=REALNATOMS.AND.r**2<MAX(a,b,c)) THEN
                                            
                            CALL ECF(.FALSE.,ECFvalue,x1,x2,y1,y2,z1,z2,px1,px2,py1,py2,pz1,pz2,a,b,c)
!                        WRITE(*,*) ECFvalue                            
                            IF(ECFvalue<1.0D0)THEN
                                            OVERLAPT2=.TRUE.
                    !                            WRITE(*,*) "overlapping ellipsoids",ECFvalue,r
                    !                                WRITE(*,'(A,3F18.10)') "Old x,r,sigma:",COORDS(K2-2,NP),r,sigma
                    !                                WRITE(*,'(A,I2,I2)') "J1,k:",J1,k
                            
                                                COORDS(K2-2,NP)=COORDS(K2-2,NP)+2*(DPRAND()-0.5D0)*2.0
                                                   COORDS(K2-1,NP)=COORDS(K2-1,NP)+2*(DPRAND()-0.5D0)*2.0
                                                COORDS(K2,NP)=COORDS(K2,NP) +2*(DPRAND()-0.5D0)*2.0
                    !                            COORDS(3*REALNATOMS+K2-2,NP)=COORDS(3*REALNATOMS+K2-2,NP)+(DPRAND()-0.5D0)*2.0*ATAN(1.0)*2
                        !                            COORDS(3*REALNATOMS+K2-1,NP)=COORDS(3*REALNATOMS+K2-1,NP)+(DPRAND()-0.5D0)*2.0*ATAN(1.0)*2
                                            
                                            

                    !                            WRITE(*,'(A)') "moving orientation"
                                                                            
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
                    
!WRITE(*,*) i,"takestep moves"


END SUBROUTINE TAKESTEPPARAM


SUBROUTINE TAKESTEPGB(NP)
! take step routine for Gay-Berne ellipsoids


USE commons
IMPLICIT NONE

INTEGER                        :: NP, J1, J2, JMAX, J3,K2,K3, JMAX2,j,k,i, REALNATOMS,JP
LOGICAL                        :: OVERLAPT, OVERLAPT2, GTEST

DOUBLE PRECISION        ::x1, y1, z1, x2, y2, z2, px1,px2,py1,py2,pz1,pz2, alphaR, &
&                       betaR,DIST(3*NATOMS/2),XMASS,YMASS,ZMASS,DMAX,VMAX,VMAX2,&
&                        S_alpha1,S_alpha2,S_alphaR,S_beta12,S_beta1R,S_beta2R,C_alpha1,C_alpha2,&
&                       C_alphaR,C_beta12,C_beta1R,C_beta2R,&
&                        VMIN,CMMAX,CMDIST(NATOMS/2),CDIST(NATOMS/2),LOCALSTEP,theta,phi,psi,PI,ECFvalue,a,b,c

DOUBLE PRECISION        :: r, DPRAND,RANDOM,&
&                         khi, khi1, DUMMY, DUMMY1, DUMMY2, DUMMY3,DUMMY4, DUMMY5, RANDOMZ


a=(GBANISOTROPYR/2.0D0)**2
b=(0.5D0)**2
c=(0.5D0)**2

!WRITE(*,*) "GBANISOTROPYR",a

PI=ATAN(1.0D0)*4
OVERLAPT=.TRUE.

REALNATOMS=NATOMS/2



    DO J1=1,REALNATOMS
        J2=3*J1
        DUMMY2=COORDS(J2-2,NP)**2+COORDS(J2-1,NP)**2+COORDS(J2,NP)**2
        IF (DUMMY2.GT.RADIUS) THEN
            WRITE(*,'(A,I5,5F20.10)') 'J1,RAD,R**2,x,y,z:', J1, RADIUS, DUMMY2, COORDS(J2-2,NP), COORDS(J2-1,NP), COORDS(J2,NP)
            PRINT*, 'initial coordinate outside container -- increase container radius'
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
!  Find the most weakly bound atom, JMAX, the second most weakly bound atom, JMAX2,
!  and the pair energy of the most tightly bound atom, VMIN. An angular step is
!  taken for JMAX if its pair energy is > ASTEP*VMIN putting the atom at a radius of
!  DMAX (or CMMAX from CM of the cluster).
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







i=0
OVERLAPT=.FALSE.
DO
i=i+1
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



                        IF(VAT(J1,NP).GT.(ASTEP(NP)*VMIN).AND.(J1.EQ.JMAX).AND.(.NOT.NORESET).AND.(J1.LT.NATOMS/2).AND.i==1) THEN
!                                    WRITE(*,*) "Angular move"
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
                        ! end of angular move block
                        ELSE IF(i==1) THEN
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
                        x1=COORDS(J2-2,NP)
                        y1=COORDS(J2-1,NP)
                        z1=COORDS(J2,NP)
                        px1=COORDS(3*REALNATOMS+J2-2,NP)
                        py1=COORDS(3*REALNATOMS+J2-1,NP)
                    pz1=COORDS(3*REALNATOMS+J2  ,NP)
                END IF
        
                DO k=J1+1,NATOMS
                    K2=3*k
                
                    IF(k<=REALNATOMS) THEN    
                            x2=COORDS(K2-2,NP)
                        y2=COORDS(K2-1,NP)
                        z2=COORDS(K2,NP)
                        px2=COORDS(3*REALNATOMS+K2-2,NP)
                        py2=COORDS(3*REALNATOMS+K2-1,NP)
                        pz2=COORDS(3*REALNATOMS+K2  ,NP)
                            r=DSQRT((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
                    END IF
                    
                    
                                    
                    IF(k<=REALNATOMS.AND.J1<=REALNATOMS.AND.r<=MAX(GBANISOTROPYR,1.0D0)) THEN
                        
                            
                            CALL ECF(.FALSE.,ECFvalue,x1,x2,y1,y2,z1,z2,px1,px2,py1,py2,pz1,pz2,a,b,c)
                                        
                            IF(ECFvalue<1.0D0)THEN
                                            OVERLAPT2=.TRUE.
                    !                            WRITE(*,*) "overlapping ellipsoids",ECFvalue,r
                    !                                WRITE(*,'(A,3F18.10)') "Old x,r,sigma:",COORDS(K2-2,NP),r,sigma
                    !                                WRITE(*,'(A,I2,I2)') "J1,k:",J1,k
                            
                                                COORDS(K2-2,NP)=COORDS(K2-2,NP)+2*(DPRAND()-0.5D0)*2.0
                                                   COORDS(K2-1,NP)=COORDS(K2-1,NP)+2*(DPRAND()-0.5D0)*2.0
                                                COORDS(K2,NP)=COORDS(K2,NP) +2*(DPRAND()-0.5D0)*2.0
                    !                            COORDS(3*REALNATOMS+K2-2,NP)=COORDS(3*REALNATOMS+K2-2,NP)+(DPRAND()-0.5D0)*2.0*ATAN(1.0)*2
                        !                            COORDS(3*REALNATOMS+K2-1,NP)=COORDS(3*REALNATOMS+K2-1,NP)+(DPRAND()-0.5D0)*2.0*ATAN(1.0)*2
                                            
                                            

                    !                            WRITE(*,'(A)') "moving orientation"
                                                                            
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
                    
!WRITE(*,*) i



END SUBROUTINE TAKESTEPGB

!
!	Dwaipayan's stuff in here
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
      DOUBLE PRECISION CTI, STI, CPI, SPI, GI(3), GJ(3), t1, p1
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
            t1            = THETA(I)
            p1            = PHI(I)
            X(OFFSET+J-2) = DCOS(p1)*DSIN(t1)
            X(OFFSET+J-1) = DSIN(p1)*DSIN(t1)
            X(OFFSET+J)   = DCOS(t1)

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

!     CALCULATE $\alpha$, $\beta$ AND $\gamma$

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

!     CALCULATE $\epsilon$

               EPS1   = DSQRT(FCT1*FCT2)
               EPS2   = 1.D0-0.5D0*GBCHIPRM*(APB*FCT7+AMB*FCT8)
               EPS    = EPSNOT*EPS1**GBNU*EPS2**GBMU

!     CALCULATE $(\sigma/\sigma_{0})^3$

               SCSIG  = 1.d0/DSQRT(1.D0-0.5D0*GBCHI*(APB*FCT3+AMB*FCT4))
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


      SUBROUTINE AAtoEuler(px,py,pz,phi,theta,chi)

      USE COMMONS

      IMPLICIT NONE
      
      INTEGER          :: I, J, K
      DOUBLE PRECISION :: R(3), P(3), A(3), RM(3,3),PI
      DOUBLE PRECISION, intent(out) :: PHI, THETA, CHI
      DOUBLE PRECISION, intent(in) :: px, py, pz
 
      PI = 4.D0 * DATAN(1.D0)

            P(1) = px
            P(2) = py
            P(3) = pz

            CALL ELLIPSOIDROTATION(P, RM)

            PHI   = DATAN2(RM(2,3),RM(1,3)) 
            IF (PHI <= 0.D0) PHI = PHI + 2.D0*PI
            
            THETA = DACOS(RM(3,3))

            CHI   = - DATAN2(RM(3,2),RM(3,1))     
            IF (CHI <= 0.D0) CHI = CHI + 2.D0*PI

            PHI    = PHI*180.D0/PI
            THETA  = THETA*180.D0/PI
            CHI    = CHI*180.D0/PI
           
!            WRITE(3,'(a5,2x,3f20.10,2x,a8,6f20.10)')                 & 
!                 'O', R(1), R(2), R(3),                              &
!                 'ellipse', 2.D0*A(1), 2.D0*A(2), 2.D0*A(3), PHI, THETA, CHI 


      END SUBROUTINE AAtoEuler

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE ELLIPSOIDROTATION (P, ROTMAT)
! gives back the three Euler angles needed for displaying an ellipsoid in Xmakemol
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
      use pymodule, only : angle,angle2,twopi

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
! sf344 additions
      LOGICAL          :: DISSOCIATED(NATOMS/2),DISSOC
      INTEGER          :: CLOSESTATOMINDEX, K1, K2, k
      DOUBLE PRECISION :: MINDISTANCE, x1, x2, y1, y2, z1, z2, RCUTARRAY(NATOMS/2), RVEC(NATOMS/2,NATOMS/2)
!      DOUBLE PRECISION, ALLOCATABLE :: RVEC(:)

!      WRITE(*,*) 'NATOMS in takestepelpsd ', NATOMS
!      IF(ALLOCATED(RVEC)) DEALLOCATE(RVEC)

      PI         = 4.D0*ATAN(1.0D0)
      IF (GBT .OR. GBDT) THEN

         RCUT       = GBEPSNOT*MAX(GBKAPPA,1.D0)

      ELSE IF (PYGPERIODICT.OR.PYBINARYT) THEN
        DO J1=1,NATOMS/2
         RCUTARRAY(J1) = 2.0D0 * MAXVAL(PYA1bin(J1,:)) 
        END DO
         RCUT = MAXVAL(RCUTARRAY)
      ELSE IF (LJCAPSIDT) THEN
        PYA1bin(:,:)=0.5D0
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
!ensure x component of particle 1 vector is within BoxLx/2 of zero. If it isn't then subtract integer number of boxlx's such that it is.
                COORDS(J2-2,NP)=COORDS(J2-2,NP)-BOXLX*NINT(COORDS(J2-2,NP)/BOXLX)
        ENDIF

        IF (PARAMONOVPBCY) THEN
!ensure y component of particle 1 vector is within BoxLy/2 of zero. If it isn't then subtract integer number of boxly's such that it is.
                COORDS(J2-1,NP)=COORDS(J2-1,NP)-BOXLY*NINT(COORDS(J2-1,NP)/BOXLY)
        END IF

        IF (PARAMONOVPBCZ) THEN
!ensure z component of particle 1 vector is within BoxLz/2 of zero. If it isn't then subtract integer number of boxlz's such that it is.
                COORDS(J2,NP)=COORDS(J2,NP)-BOXLZ*NINT(COORDS(J2,NP)/BOXLZ)
        ENDIF

         DUMMY2 = COORDS(J2-2,NP)**2 + COORDS(J2-1,NP)**2 + COORDS(J2,NP)**2

         IF (DUMMY2 .GT. RADIUS) THEN
! bring back the molecule within the radius
                COORDS(J2-2,NP)=COORDS(J2-2,NP)-SQRT(RADIUS)*NINT(COORDS(J2-2,NP)/SQRT(RADIUS))
                COORDS(J2-1,NP)=COORDS(J2-1,NP)-SQRT(RADIUS)*NINT(COORDS(J2-1,NP)/SQRT(RADIUS))
                COORDS(J2,NP)=COORDS(J2,NP)-SQRT(RADIUS)*NINT(COORDS(J2,NP)/SQRT(RADIUS))
    WRITE(MYUNIT,'(A,2F20.10)') 'initial coordinate outside container -- bringing molecule back within the container radius'

!            WRITE(*,'(A,I5,5F20.10)') 'J1,RAD,R**2,x,y,z:', J1, RADIUS, DUMMY2, COORDS(J2-2,NP), &
!                                       COORDS(J2-1,NP), COORDS(J2,NP)
!            PRINT*, 'initial coordinate outside container -- increase container radius'
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

!     Find the most weakly bound atom, JMAX, the second most weakly bound atom, JMAX2,
!     and the pair energy of the most tightly bound atom, VMIN. An angular step is
!     taken for JMAX if its pair energy is > ASTEP*VMIN putting the atom at a radius of
!     DMAX (or CMMAX from CM of the cluster).
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
!        WRITE(*,*) 'in here'



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

            CALL ELPSDRTMT (P1, PYA1bin(J1,:), AE)

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

               CALL ELPSDRTMT (P2, PYA1bin(J2,:), BE)

               RIJ    = RI - RJ
               ABSRIJ = DSQRT(DOT_PRODUCT(RIJ,RIJ))
                        
               IF (ABSRIJ < RCUT) THEN 

!     DETERMINE ELLIPTIC CONTACT FUNCTION

                  CALL BRENTMIN (0.D0, 0.51D0, 1.D0, AE, BE, RIJ, XMIN, FMIN)

                  ECFVAL = - FMIN
! allow for a slight overlap 
                  IF (ECFVAL < PYOVERLAPTHRESH) THEN
!                        WRITE(*,*) 'atoms overlapping', J1, J2, ECFVAL, ABSRIJ
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


! now determine whether any particle has dissociated
 IF(PYGPERIODICT) THEN
  DO
  DISSOC=.FALSE.
   DO J1=1,natoms/2-1
    J2=3*J1
    x1=COORDS(J2-2,NP)
    y1=COORDS(J2-1,NP)
    z1=COORDS(J2,NP)
   
    DISSOCIATED(J1)=.TRUE.
    DO K1=1,natoms/2
      K2=3*K1
      x2=COORDS(K2-2,NP)
      y2=COORDS(K2-1,NP)
      z2=COORDS(K2,NP)
      RVEC(J1,K1)=SQRT((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
!      WRITE(*,*) J1,K1,RVEC(J1,K1), 4*RCUT
      IF(RVEC(J1,K1)<4*RCUT) DISSOCIATED(J1)=.FALSE.
!      WRITE(*,*) 'in here, cutoff=',pcutoff,J1
    END DO
    IF(DISSOCIATED(J1)) THEN
        ! determine the closest neighbor
        MINDISTANCE = HUGE(1.0D0)
        DO K1=1,natoms/2
         IF(K1==J1) CYCLE
         IF(RVEC(J1,K1)<MINDISTANCE) THEN
           MINDISTANCE=RVEC(J1,K1)
           CLOSESTATOMINDEX=K1
         END IF
        END DO
        WRITE(MYUNIT,'(A,I6,A,F10.6,I6)') 'Atom ', J1, ' is dissociated, distance to the closest atom: ', &
        & MINDISTANCE, CLOSESTATOMINDEX
        DISSOC=.TRUE.
        ! move the dissociated atom closer to the closest atom
        COORDS(J2-2,NP)=0.5D0*(COORDS(3*CLOSESTATOMINDEX-2,NP)+COORDS(J2-2,NP))
        COORDS(J2-1,NP)=0.5D0*(COORDS(3*CLOSESTATOMINDEX-1,NP)+COORDS(J2-1,NP))
        COORDS(J2  ,NP)=0.5D0*(COORDS(3*CLOSESTATOMINDEX  ,NP)+COORDS(J2  ,NP))
    END IF
   END DO
   IF(.NOT.DISSOC) EXIT
  END DO
        DO J1=NATOMS/2+1,NATOMS
            J2=3*J1
            angle2=dot_product(COORDS(J2-2:J2,NP),COORDS(J2-2:J2,NP))

            if(angle2>twopi**2) then
                angle2=sqrt(angle2)
                angle = angle2 - int(angle2/twopi)*twopi
                COORDS(J2-2:J2,NP)=COORDS(J2-2:J2,NP)/angle2 * angle
            end if
        END DO

RETURN
!     CHECK FOR OVERLAP AGAIN (not working!)
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
use commons, only : NATOMS, COORDS, PYBINARYTYPE1,MYUNIT,SWAPMOVEST,PYSWAP

implicit none

integer :: NP, RANDOM1, RANDOM2, J1, J2
double precision :: COORDSSTORE(3,4),DPRAND

! get two random integer numbers, one for atom index type 1, one for atom index type 2

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
        ! select a body from the two types
        DO J1=1,3
         COORDSSTORE(J1,1)=COORDS(3*(RANDOM1-1)+J1,NP)
         COORDSSTORE(J1,2)=COORDS(3*(RANDOM2-1)+J1,NP)
         COORDSSTORE(J1,3)=COORDS(3*NATOMS/2+3*(RANDOM1-1)+J1,NP)
         COORDSSTORE(J1,4)=COORDS(3*NATOMS/2+3*(RANDOM2-1)+J1,NP)
        END DO

        ! swap coordinates
        WRITE(MYUNIT,*) 'Swapping atoms ', RANDOM1, RANDOM2
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
              WRITE(MYUNIT,*) 'all atoms swapped, restarting'
!      SWAPMOVEST=.FALSE.  
              PYSWAP(1)=1
              PYSWAP(2)=PYBINARYTYPE1+1
              EXIT
        END IF
END DO
END SUBROUTINE TAKESTEPSWAPMOVES

SUBROUTINE PYPES(XORIG)
use commons,only : natoms
implicit none

integer :: j0,j1,j2,j3,nsteps
double precision :: interval, XORIG(3*NATOMS),XNEW(3*NATOMS),GRAD(3*NATOMS)
double precision, allocatable   :: PES(:,:,:),X1(:),Y1(:),Z1(:)

interval=0.2D0
nsteps=60

XNEW(:)=XORIG(:)

ALLOCATE(PES(nsteps+1,nsteps+1,nsteps+1),X1(nsteps+1),Y1(nsteps+1),Z1(nsteps+1))
open(unit=955,file='py_pes',status='unknown')
do j0=0,nsteps
 XNEW(2)=XORIG(2)+j0*interval
 Y1(j0+1)=XNEW(2)
 do j1=0,nsteps
   WRITE(*,*) 'working... ', j0,j1
   XNEW(1)=XORIG(1)+j1*interval
   X1(j1+1)=XNEW(1)
   do j2=0,nsteps
      XNEW(3)=XORIG(3)+j2*interval
      Z1(j2+1)=XNEW(3)
      CALL PYGPERIODIC(XNEW,GRAD,PES(j0+1,j1+1,j2+1),.TRUE.)
      IF(PES(j0+1,j1+1,j2+1)>1.0D9) THEN
        PES=1.0D9
      ELSE IF(PES(j0+1,j1+1,j2+1)<-1.0D3) THEN
        PES=-1.0D9
      END IF
         write(955,'(4F20.5)') X1(j1+1),Y1(j0+1),Z1(j2+1),PES(j0+1,j1+1,j2+1)
   end do
 end do
end do

do j0=1,nsteps+1
 do j1=1,nsteps+1
   do j2=1,nsteps+1
!        write(955,'(4F20.5)') X1(j1),Y1(j0),Z1(j2),PES(j0,j1,j2)
   end do
 end do
end do
close(955)
END SUBROUTINE PYPES

SUBROUTINE INITIALISELJCAPSIDMODEL

use commons, only: BOXLX,BOXLY,BOXLZ,PARAMONOVPBCX,PARAMONOVPBCY,PARAMONOVPBCZ,PCUTOFF,PARAMONOVCUTOFF,&
                &       natoms,pya1bin,pya2bin,pysignot,pyepsnot,radift,LJSITE,BLJSITE,PEPSILON1,&
                &       PSCALEFAC1,PSCALEFAC2,MAXINTERACTIONS,PYBINARYT,PYBINARYTYPE1,MYUNIT,VT
use ljcapsidmodule 

implicit none
   WRITE(MYUNIT,*) 'initialising variables for LJ capsid model',NATOMS 
! allocate arrays

    ALLOCATE(RMIvec(natoms/2,3,3),DPI1RMvec(natoms/2,3,3), DPI2RMvec(natoms/2,3,3), DPI3RMvec(natoms/2,3,3))
    ALLOCATE(PSCALEFAC1vec(natoms/2),PSCALEFAC2vec(natoms/2))
    ALLOCATE(AEZR1(NATOMS/2,3,3), AEZR2(NATOMS/2,3,3))
    ALLOCATE(epsilon1(4,natoms/2,natoms/2))
      I3(:,:)    = 0.D0
      AEZR1(:,:,:) = 0.D0
      AEZR2(:,:,:) = 0.D0



     pi=ATAN(1.0D0)*4.0
     twopi=2.0D0*pi

       epsilon1(1:3,:,:)=PEPSILON1(1)   ! we are going to use epsilon1 for the extra LJ sites
       epsilon1(4,:,:)=PYEPSNOT
       sigma1=1.0D0 ! PSIGMA1 is nonexistent from now on
       PSCALEFAC1vec(:)=PSCALEFAC1(1)
       PSCALEFAC2vec(:)=PSCALEFAC2(1)

      REALNATOMS = NATOMS/2
      OFFSET     = 3*REALNATOMS
      
      DO K1 = 1, 3
        I3(K1,K1) = 1.0D0
      ENDDO
END SUBROUTINE INITIALISELJCAPSIDMODEL

SUBROUTINE LJCAPSIDMODEL (X, G, ENERGY, GTEST)
use commons, only:  natoms,pysignot,pyepsnot,LJSITE,BLJSITE,PEPSILON1,&
                &       PSCALEFAC1,PSCALEFAC2,MAXINTERACTIONS,PYBINARYT,PYBINARYTYPE1,MYUNIT,VT
use ljcapsidmodule

implicit none
double precision :: dummy, vdummy, dvdummy(12), dummy1, dummy2, p(3), term2(maxinteractions), term3(maxinteractions), &
                  & xlj(maxinteractions,2,3), rljvec(maxinteractions,3), rljunitvec(maxinteractions,3), vecsbf(3),&
                  & drlj(maxinteractions,12), rlj(maxinteractions), rlj2(maxinteractions), &
                  & llj(12,maxinteractions), dllj1(maxinteractions,12),attr(maxinteractions), &
                  & tij(3), tji(3), fij(3), energy, x(3*natoms), g(3*natoms), sigma1vec(maxinteractions), &
                  & RMI(3,3), RMJ(3,3), RI(3), RJ(3), &
                  & DPI1RM(3,3), DPI2RM(3,3), DPI3RM(3,3), DPJ1RM(3,3), DPJ2RM(3,3), DPJ3RM(3,3)
integer          :: k
!integer          :: j1, j3, j5, k, k1, realnatoms
logical          :: gtest


VT(1:NATOMS/2)=0.0D0
 
term2(:)=1.0D0
term3(:)=0.0D0
attr(1:3)=0.0D0
attr(4)=1.0D0
sigma1vec(1:3)=1.0
sigma1vec(4)=1.0
!pi=ATAN(1.0D0)*4.0
!twopi=2.0D0*pi

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
            angle=sqrt(dot_product(P,P))
            if(angle>twopi) then
! normalise angle-axis coordinates
                X(J5-2:J5)=X(J5-2:J5)/angle
                do
                  angle=angle-twopi
                  if(angle<2*pi) exit
                end do
! multiply with new angle
                X(J5-2:J5)=X(J5-2:J5)*angle
            end if
!        WRITE(*,*) 'calling RMDRVT', P, RMIvec(J1,:,:), DPI1RMvec(J1,:,:), DPI2RMvec(J1,:,:), DPI3RMvec(J1,:,:), GTEST
            CALL RMDRVT(P, RMIvec(J1,:,:), DPI1RMvec(J1,:,:), DPI2RMvec(J1,:,:), DPI3RMvec(J1,:,:), GTEST)

        END DO

         DO J1 = 1, REALNATOMS - 1

            J3      = 3*J1
            J5      = OFFSET + J3
            RI      = X(J3-2:J3)
            P       = X(J5-2:J5)
!     ROTATION MATRIX

!            CALL RMDRVT(P, RMI, DPI1RM, DPI2RM, DPI3RM, GTEST)
            RMI(:,:)=RMIvec(J1,:,:)
            DPI1RM(:,:)=DPI1RMvec(J1,:,:)
            DPI2RM(:,:)=DPI2RMvec(J1,:,:)
            DPI3RM(:,:)=DPI3RMvec(J1,:,:)


           DO J2 = J1 + 1, REALNATOMS
               J4     = 3*J2
               J6     = OFFSET + J4
               RJ     = X(J4-2:J4)
               P      = X(J6-2:J6)

               RMJ(:,:)=RMIvec(J2,:,:)
               DPJ1RM(:,:)=DPI1RMvec(J2,:,:)
               DPJ2RM(:,:)=DPI2RMvec(J2,:,:)
               DPJ3RM(:,:)=DPI3RMvec(J2,:,:)

               VDUMMY=0.0D0
               dVDUMMY(:)=0.0D0
            do k=1,maxinteractions ! k=1 -- interaction between repulsive primary 'apex' sites
                         ! k=2 and k=3 -- interaction between secondary and primary 'apex' sites 
                         ! k=4  -- LJ site at the centre (normal 12-6 interaction)
                        vecsbf(1)=1.0D0
                        vecsbf(2)=0.0D0
                        vecsbf(3)=0.0D0
                        ! trying to modify code to allow for binary systems. 
                        ! apex site heights will be defined in absolute units,
                        ! hence PYA1bin(J1,1) etc. will be removed from below

                   IF(k==1) THEN
                        DUMMY1=PSCALEFAC1vec(J1)!*PYA1bin(J1,1)
                        DUMMY2=PSCALEFAC1vec(J2)!*PYA1bin(J2,1)
                   ELSE IF(k==2) THEN
                        DUMMY1=PSCALEFAC1vec(J1)!*PYA1bin(J1,1)
                        DUMMY2=-PSCALEFAC2vec(J2)!*PYA1bin(J2,1)
                   ELSE IF(k==3) THEN
                        DUMMY1=-PSCALEFAC2vec(J1)!*PYA1bin(J1,1)
                        DUMMY2=PSCALEFAC1vec(J2)!*PYA1bin(J2,1)
                   ELSE IF(k==4) THEN
                        DUMMY1=0.0D0
                        DUMMY2=0.0D0
                   END IF
                        ! first particle
                        xlj(k,1,:)=RI+DUMMY1*MATMUL(RMI,vecsbf)    ! vecsbf: (1,0,0) in the body frame of ellipsoid

                        ! second particle
                        xlj(k,2,:)=RJ+DUMMY2*MATMUL(RMJ,vecsbf)

                        ! separation between the LJ sites
                        rlj2(k)=(xlj(k,2,1)-xlj(k,1,1))**2+(xlj(k,2,2)-xlj(k,1,2))**2+(xlj(k,2,3)-xlj(k,1,3))**2
                        rlj(k)=sqrt(rlj2(k))
                        rljvec(k,1)=xlj(k,2,1)-xlj(k,1,1)
                        rljvec(k,2)=xlj(k,2,2)-xlj(k,1,2)
                        rljvec(k,3)=xlj(k,2,3)-xlj(k,1,3)

                        DUMMY=1.0D0/rlj(k)
                        rljunitvec(k,:)=rljvec(k,:)*DUMMY !/rlj(k)

                        drlj(k,1)=DUMMY*(xlj(k,2,1)-xlj(k,1,1))         !drlj/dx1
                        drlj(k,2)=DUMMY*(xlj(k,2,2)-xlj(k,1,2))         !drlj/dy1
                        drlj(k,3)=DUMMY*(xlj(k,2,3)-xlj(k,1,3))         !drlj/dz1
                        drlj(k,4)=-drlj(k,1)                               !drlj/dx2
                        drlj(k,5)=-drlj(k,2)                               !drlj/dy2
                        drlj(k,6)=-drlj(k,3)                               !drlj/dz2
                        drlj(k,7) =-DUMMY*DUMMY1*DOT_PRODUCT(rljvec(k,:),MATMUL(DPI1RM,vecsbf)) !drlj/dpx1
                        drlj(k,8) =-DUMMY*DUMMY1*DOT_PRODUCT(rljvec(k,:),MATMUL(DPI2RM,vecsbf)) !drlj/dpy1
                        drlj(k,9) =-DUMMY*DUMMY1*DOT_PRODUCT(rljvec(k,:),MATMUL(DPI3RM,vecsbf)) !drlj/dpz1
                        drlj(k,10) =  DUMMY*DUMMY2*DOT_PRODUCT(rljvec(k,:),MATMUL(DPJ1RM,vecsbf)) !drlj/dpx2
                        drlj(k,11) =  DUMMY*DUMMY2*DOT_PRODUCT(rljvec(k,:),MATMUL(DPJ2RM,vecsbf)) !drlj/dpy2
                        drlj(k,12) =  DUMMY*DUMMY2*DOT_PRODUCT(rljvec(k,:),MATMUL(DPJ3RM,vecsbf)) !drlj/dpz2

              ! interaction between the LJ sites:
                        
                        LLJ(1,k)=sigma1vec(k)*DUMMY !/rlj(k)
                        LLJ(2,k)=LLJ(1,k)**2
                        LLJ(3,k)=LLJ(2,k)*LLJ(1,k)
                        LLJ(4,k)=LLJ(2,k)**2
                        LLJ(5,k)=LLJ(4,k)*LLJ(1,k)
                        LLJ(6,k)=LLJ(4,k)*LLJ(2,k)
                        LLJ(7,k)=LLJ(6,k)*LLJ(1,k)
                        LLJ(11,k)=LLJ(5,k)*LLJ(6,k)
                        LLJ(12,k)=LLJ(6,k)*LLJ(6,k)

!                            DUMMY=1.0D0/rlj(k)
!                            DUMMY=DUMMY**2
                            DO j=1,12
                                dLLJ1(k,j) =-sigma1vec(k)*DUMMY*DUMMY*drlj(k,j)
                            END DO
                 VDUMMY=VDUMMY+4.0D0*epsilon1(k,J1,J2)*term2(k)*(LLJ(12,k) - attr(k)*LLJ(6,k))
                                     ! extra LJ sites, repulsive, and LJ site in the middle (attr(4)=1, attr(1:3)=0)
!                 IF(k==3) WRITE(*,*) 'repulsive energies ', 4.0D0*epsilon1(k,J1,J2)*LLJ(12,k)*term2(k)
             do j=1,3
               dVDUMMY(j) = dVDUMMY(j) + 4.0D0*epsilon1(k,J1,J2)*(12.0D0*LLJ(11,k)*dLLJ1(k,j))*term2(k) + &
                          & 4.0D0*epsilon1(k,J1,J2)*LLJ(12,k)*term3(k)*drlj(k,j)! extra LJ site derivatives, currently only repulsive
               dVDUMMY(j) = dVDUMMY(j) - attr(k)*(4.0D0*epsilon1(k,J1,J2)*(6.0D0*LLJ(5,k)*dLLJ1(k,j))*term2(k) + &
                          & 4.0D0*epsilon1(k,J1,J2)*LLJ(6,k)*term3(k)*drlj(k,j))
               FIJ(j) = dVDUMMY(j)
             end do
             do j=7,12
               dVDUMMY(j) = dVDUMMY(j) + 4.0D0*epsilon1(k,J1,J2)*(12.0D0*LLJ(11,k)*dLLJ1(k,j))*term2(k) + &
                          & 4.0D0*epsilon1(k,J1,J2)*LLJ(12,k)*term3(k)*drlj(k,j)! extra LJ site derivatives, currently only repulsive
               dVDUMMY(j) = dVDUMMY(j) - attr(k)*(4.0D0*epsilon1(k,J1,J2)*(6.0D0*LLJ(5,k)*dLLJ1(k,j))*term2(k) + &
                          & 4.0D0*epsilon1(k,J1,J2)*LLJ(6,k)*term3(k)*drlj(k,j))
             end do
            end do !k=1,4

               TIJ = dVDUMMY(7:9)
               TJI = dVDUMMY(10:12)

!                G(J3-2)=G(J3-2)+dVDUMMY(1)
!                G(j4-2)=G(j4-2)+dVDUMMY(4)
!                G(J3-1)=G(J3-1)+dVDUMMY(2)
!                G(j4-1)=G(j4-1)+dVDUMMY(5)
!                G(J3  )=G(J3  )+dVDUMMY(3)
!                G(j4  )=G(j4  )+dVDUMMY(6)



!              G(J3-2:J3) = G(J3-2:J3) + dvdummy(1:3)
!              G(J4-2:J4) = G(J4-2:J4) + dvdummy(4:6)
!              G(J5-2:J5) = G(J5-2:J5) + dvdummy(7:9)
!              G(J6-2:J6) = G(J6-2:J6) + dvdummy(10:12)

               ENERGY = ENERGY + VDUMMY
               VT(J1) = VT(J1) + VDUMMY
               VT(J2) = VT(J2) + VDUMMY        ! pair potentials

               G(J3-2:J3) = G(J3-2:J3) - FIJ
               G(J4-2:J4) = G(J4-2:J4) + FIJ
               G(J5-2:J5) = G(J5-2:J5) + TIJ
               G(J6-2:J6) = G(J6-2:J6) + TJI

        end do
       end do

END SUBROUTINE LJCAPSIDMODEL

