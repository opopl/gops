      SUBROUTINE MULTISITEPY (X, G, ENERGY, GTEST)

      USE COMMONS, ONLY: NATOMS, NPYSITE, PYA1, PYA2, PYSIGNOT, PYEPSNOT, RADIFT 

      IMPLICIT NONE

      INTEGER          :: I, J, J1, J2, J3, J4, J5, J6, K1, K2, REALNATOMS, OFFSET
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS)
      DOUBLE PRECISION :: PST(NATOMS,3), EZRI1(3,3), EZRI2(3,3), EZRJ1(3,3), EZRJ2(3,3)
      DOUBLE PRECISION :: RI(3), RJ(3), RIJ(3), NR(3), RIJSQ, ABSRIJ, P(3), THETA, THETA2, CT, ST
      DOUBLE PRECISION :: I3(3,3), AE1(3,3), BE1(3,3), AE2(3,3), BE2(3,3), APB(3,3), APBINV(3,3)
      DOUBLE PRECISION :: RSTI(3), RSTJ(3), FCNT1, FCNT2, SRTFI1, SRTFI2, FMIN, LAMDAC, ENERGY
      DOUBLE PRECISION :: RHO1, RHO1SQ, RHO16, RHO112, RHO2, RHO2SQ, RHO26
      DOUBLE PRECISION :: FCTR1, FCTR2, DV1DF1, DV2DF2, DV1DR, DV2DR
      DOUBLE PRECISION :: DRDPI1, DRDPI2, DRDPI3, DRDPJ1, DRDPJ2, DRDPJ3 
      DOUBLE PRECISION :: DF1PI1, DF1PI2, DF1PI3, DF2PI1, DF2PI2, DF2PI3
      DOUBLE PRECISION :: DF1PJ1, DF1PJ2, DF1PJ3, DF2PJ1, DF2PJ2, DF2PJ3 
      DOUBLE PRECISION :: RMI(3,3), RMJ(3,3), E(3,3)
      DOUBLE PRECISION :: DRMI1(3,3), DRMI2(3,3), DRMI3(3,3), DRMJ1(3,3), DRMJ2(3,3), DRMJ3(3,3)
      DOUBLE PRECISION :: DF1DR(3), DF2DR(3), DG1DR(3), DG2DR(3)
      DOUBLE PRECISION :: ARIBRJ(3), XC(3), XCMRI(3), XCMRJ(3), FIJ(3), TIJ(3), TJI(3)
      DOUBLE PRECISION :: D1ABEZ(3,3), D2ABEZ(3,3), D3ABEZ(3,3), D1ABE(3,3), D2ABE(3,3), D3ABE(3,3) 
      LOGICAL          :: GTEST

      I3(:,:)      = 0.D0
       NPYSITE = 100
!     FROM INPUT PARAMETERS

      REALNATOMS = NATOMS/2
      OFFSET     = 3*REALNATOMS

      DO K1 = 1, 3

         I3(K1,K1)      = 1.D0

      ENDDO

      CALL DEFMSPY(PST)

      IF (GTEST) THEN

         ENERGY = 0.D0
         G(:)   = 0.D0

         DO J1 = 1, REALNATOMS - 1

            J3      = 3*J1
            J5      = OFFSET + J3
            RI      = X(J3-2:J3)
            P       = X(J5-2:J5)

!     ROTATION MATRIX

            CALL RMDRVT(P, RMI, DRMI1, DRMI2, DRMI3, GTEST)

!     BEGIN INNER LOOP OVER PARTICLES

            DO J2 = J1 + 1, REALNATOMS

               J4     = 3*J2
               J6     = OFFSET + J4
               RJ     = X(J4-2:J4) 
               P      = X(J6-2:J6)

!     ROTATION MATRIX

               CALL RMDRVT(P, RMJ, DRMJ1, DRMJ2, DRMJ3, GTEST)

               DO I = 1, NPYSITE

!     OBTAIN THE SITE POSITION RELATIVE TO THE COM IN THE SPACE-FIXED FRAME

                  RSTI(:) = MATMUL(RMI,PST(I,:))
                  
                  CALL SITEBF (I, EZRI1, EZRI2)

                  AE1 = MATMUL(RMI,(MATMUL(EZRI1(:,:),(TRANSPOSE(RMI)))))

                  IF (RADIFT) THEN

                     AE2 = MATMUL(RMI,(MATMUL(EZRI2(:,:),(TRANSPOSE(RMI)))))

                  ENDIF
               
                  DO J = 1, NPYSITE

                     RSTJ(:) = MATMUL(RMJ,PST(J,:))

                     CALL SITEBF (J, EZRJ1, EZRJ2)

                     BE1 = MATMUL(RMJ,(MATMUL(EZRJ1(:,:),(TRANSPOSE(RMJ)))))

                     IF (RADIFT) THEN
   
                        BE2 = MATMUL(RMJ,(MATMUL(EZRJ2(:,:),(TRANSPOSE(RMJ)))))

                     ENDIF

!     CALCULATE SEPARATION

                     RIJ    = RI - RJ + RSTI - RSTJ
                     RIJSQ  = DOT_PRODUCT(RIJ,RIJ)
                     ABSRIJ = DSQRT(RIJSQ)
                     NR     = RIJ / ABSRIJ

!     CALCULATE ECF
                     CALL BRENTMIN (0.D0, 0.51D0, 1.D0, AE1, BE1, RIJ, LAMDAC, FMIN)

                     FCNT1   = - FMIN
                     SRTFI1  = 1.D0 / DSQRT(FCNT1)
                     APB     = LAMDAC * AE1 + (1.D0 - LAMDAC) * BE1

                     CALL MTRXIN (APB, APBINV)

                     ARIBRJ =  LAMDAC * MATMUL(AE1,(RI+RSTI)) + (1.D0 - LAMDAC) * MATMUL(BE1,(RJ+RSTJ))
                     XC     =  MATMUL(APBINV, ARIBRJ)
                     XCMRI  = XC - RI - RSTI
                     XCMRJ  = XC - RJ - RSTJ
                     DF1DR  = - 2.D0 * LAMDAC * MATMUL(AE1,XCMRI)

                     D1ABEZ = MATMUL(DRMI1,EZRI1(:,:))
                     D1ABE  = MATMUL(D1ABEZ,TRANSPOSE(RMI)) + MATMUL(RMI,(TRANSPOSE(D1ABEZ)))

                     D2ABEZ = MATMUL(DRMI2,EZRI1(:,:))
                     D2ABE  = MATMUL(D2ABEZ,(TRANSPOSE(RMI))) + MATMUL(RMI,(TRANSPOSE(D2ABEZ)))

                     D3ABEZ = MATMUL(DRMI3,EZRI1(:,:))
                     D3ABE  = MATMUL(D3ABEZ,(TRANSPOSE(RMI))) + MATMUL(RMI,(TRANSPOSE(D3ABEZ)))

                     DF1PI1 = LAMDAC*(DOT_PRODUCT(XCMRI,MATMUL(D1ABE,XCMRI)) &
                            - 2.D0*DOT_PRODUCT(XCMRI,MATMUL(AE1,MATMUL(DRMI1,PST(I,:))))) 
                     DF1PI2 = LAMDAC*(DOT_PRODUCT(XCMRI,MATMUL(D2ABE,XCMRI)) &
                            - 2.D0*DOT_PRODUCT(XCMRI,MATMUL(AE1,MATMUL(DRMI2,PST(I,:)))))
                     DF1PI3 = LAMDAC*(DOT_PRODUCT(XCMRI,MATMUL(D3ABE,XCMRI)) &
                            - 2.D0*DOT_PRODUCT(XCMRI,MATMUL(AE1,MATMUL(DRMI3,PST(I,:)))))

                     D1ABEZ = MATMUL(DRMJ1,EZRJ1(:,:))
                     D1ABE  = MATMUL(D1ABEZ,(TRANSPOSE(RMJ))) + MATMUL(RMJ,(TRANSPOSE(D1ABEZ)))

                     D2ABEZ = MATMUL(DRMJ2,EZRJ1(:,:))
                     D2ABE  = MATMUL(D2ABEZ,(TRANSPOSE(RMJ))) + MATMUL(RMJ,(TRANSPOSE(D2ABEZ)))

                     D3ABEZ = MATMUL(DRMJ3,EZRJ1(:,:))
                     D3ABE  = MATMUL(D3ABEZ,(TRANSPOSE(RMJ))) + MATMUL(RMJ,(TRANSPOSE(D3ABEZ))) 
               
                     DF1PJ1 = (1.D0-LAMDAC)*(DOT_PRODUCT(XCMRJ,MATMUL(D1ABE,XCMRJ)) &
                            - 2.D0*DOT_PRODUCT(XCMRJ,MATMUL(BE1,MATMUL(DRMJ1,PST(J,:)))))
                     DF1PJ2 = (1.D0-LAMDAC)*(DOT_PRODUCT(XCMRJ,MATMUL(D2ABE,XCMRJ)) &
                            - 2.D0*DOT_PRODUCT(XCMRJ,MATMUL(BE1,MATMUL(DRMJ2,PST(J,:)))))
                     DF1PJ3 = (1.D0-LAMDAC)*(DOT_PRODUCT(XCMRJ,MATMUL(D3ABE,XCMRJ)) &
                            - 2.D0*DOT_PRODUCT(XCMRJ,MATMUL(BE1,MATMUL(DRMJ3,PST(J,:))))) 

                     RHO1   = PYSIGNOT / (ABSRIJ - ABSRIJ*SRTFI1 + PYSIGNOT)
                     RHO1SQ = RHO1*RHO1
                     RHO16  = RHO1SQ*RHO1SQ*RHO1SQ
                     RHO112 = RHO16 * RHO16

                     FCTR1  = 0.5D0*ABSRIJ*SRTFI1/(FCNT1*PYSIGNOT)
                     DG1DR  = (1.D0-SRTFI1)*NR/PYSIGNOT + FCTR1*DF1DR
                     DV1DF1 = -2.D0*RHO112*RHO1*FCTR1
                     DV1DR  = -2.D0*RHO112*RHO1*(1.D0-SRTFI1)/PYSIGNOT

                     IF (RADIFT) THEN

                        CALL BRENTMIN (0.D0, 0.51D0, 1.D0, AE2, BE2, RIJ, LAMDAC, FMIN)

                        FCNT2   = - FMIN
                        SRTFI2  = 1.D0 / DSQRT(FCNT2)
                        APB     = LAMDAC * AE2 + (1.D0 - LAMDAC) * BE2

                        CALL MTRXIN (APB, APBINV)

                        ARIBRJ = LAMDAC * MATMUL(AE2,(RI+RSTI)) + (1.D0 - LAMDAC) * MATMUL(BE2,(RJ+RSTJ))
                        XC     = MATMUL(APBINV, ARIBRJ)
                        XCMRI  = XC - RI - RSTI
                        XCMRJ  = XC - RJ - RSTJ
                        DF2DR  = - 2.D0 * LAMDAC * MATMUL(AE2,XCMRI)

                        RHO2   = PYSIGNOT / (ABSRIJ - ABSRIJ*SRTFI2 + PYSIGNOT)
                        RHO2SQ = RHO2*RHO2
                        RHO26  = RHO2SQ*RHO2SQ*RHO2SQ
               
                        FCTR2  = 0.5D0*ABSRIJ*SRTFI2/(FCNT2*PYSIGNOT)
                        DG2DR  = (1.D0-SRTFI2)*NR/PYSIGNOT+FCTR2*DF2DR
                        DV2DF2 = RHO26*RHO2*FCTR2
                        DV2DR  = RHO26*RHO2*(1.D0-SRTFI2)/PYSIGNOT

                        D1ABEZ = MATMUL(DRMI1,EZRI2(:,:))
                        D1ABE  = MATMUL(D1ABEZ,(TRANSPOSE(RMI))) + MATMUL(RMI,(TRANSPOSE(D1ABEZ)))

                        D2ABEZ = MATMUL(DRMI2,EZRI2(:,:))
                        D2ABE  = MATMUL(D2ABEZ,(TRANSPOSE(RMI))) + MATMUL(RMI,(TRANSPOSE(D2ABEZ)))

                        D3ABEZ = MATMUL(DRMI3,EZRI2(:,:))
                        D3ABE  = MATMUL(D3ABEZ,(TRANSPOSE(RMI))) + MATMUL(RMI,(TRANSPOSE(D3ABEZ)))

                        DF2PI1 = LAMDAC*(DOT_PRODUCT(XCMRI,MATMUL(D1ABE,XCMRI)) &
                               - 2.D0*DOT_PRODUCT(XCMRI,MATMUL(AE2,MATMUL(DRMI1,PST(I,:)))))
                        DF2PI2 = LAMDAC*(DOT_PRODUCT(XCMRI,MATMUL(D2ABE,XCMRI)) &
                               - 2.D0*DOT_PRODUCT(XCMRI,MATMUL(AE2,MATMUL(DRMI2,PST(I,:)))))
                        DF2PI3 = LAMDAC*(DOT_PRODUCT(XCMRI,MATMUL(D3ABE,XCMRI)) &
                               - 2.D0*DOT_PRODUCT(XCMRI,MATMUL(AE2,MATMUL(DRMI3,PST(I,:)))))
                        
                        D1ABEZ = MATMUL(DRMJ1,EZRJ2(:,:))
                        D1ABE  = MATMUL(D1ABEZ,(TRANSPOSE(RMJ))) + MATMUL(RMJ,(TRANSPOSE(D1ABEZ)))

                        D2ABEZ = MATMUL(DRMJ2,EZRJ2(:,:))
                        D2ABE  = MATMUL(D2ABEZ,(TRANSPOSE(RMJ))) + MATMUL(RMJ,(TRANSPOSE(D2ABEZ)))

                        D3ABEZ = MATMUL(DRMJ3,EZRJ2(:,:))
                        D3ABE  = MATMUL(D3ABEZ,(TRANSPOSE(RMJ))) + MATMUL(RMJ,(TRANSPOSE(D3ABEZ)))

                        DF2PJ1 = (1.D0-LAMDAC)*(DOT_PRODUCT(XCMRJ,MATMUL(D1ABE,XCMRJ)) &
                               - 2.D0*DOT_PRODUCT(XCMRJ,MATMUL(BE2,MATMUL(DRMJ1,PST(J,:)))))
                        DF2PJ2 = (1.D0-LAMDAC)*(DOT_PRODUCT(XCMRJ,MATMUL(D2ABE,XCMRJ)) &
                               - 2.D0*DOT_PRODUCT(XCMRJ,MATMUL(BE2,MATMUL(DRMJ2,PST(J,:)))))
                        DF2PJ3 = (1.D0-LAMDAC)*(DOT_PRODUCT(XCMRJ,MATMUL(D3ABE,XCMRJ)) &
                               - 2.D0*DOT_PRODUCT(XCMRJ,MATMUL(BE2,MATMUL(DRMJ3,PST(J,:)))))

                     ELSE

                        RHO2   = RHO1
                        RHO26  = RHO16
                        DG2DR  = DG1DR
                        DV2DF2 = RHO26*RHO2*FCTR1
                        DV2DR  = RHO26*RHO2*(1.D0-SRTFI1)/PYSIGNOT
                        DF2PI1 = DF1PI1
                        DF2PI2 = DF1PI2
                        DF2PI3 = DF1PI3
                        DF2PJ1 = DF1PJ1
                        DF2PJ2 = DF1PJ2
                        DF2PJ3 = DF1PJ3

                     ENDIF
             
!     CALCULATE PY POTENTIAL ENERGY

!     CALCULATE GRADIENT

                     FIJ = 2.D0*RHO112*RHO1*DG1DR - RHO26*RHO2*DG2DR
                    
                     DRDPI1 = DOT_PRODUCT(NR,MATMUL(DRMI1,PST(I,:)))
                     DRDPI2 = DOT_PRODUCT(NR,MATMUL(DRMI2,PST(I,:)))
                     DRDPI3 = DOT_PRODUCT(NR,MATMUL(DRMI3,PST(I,:)))

                     TIJ(1) = DV1DF1*DF1PI1 + DV1DR*DRDPI1 + DV2DF2*DF2PI1 + DV2DR*DRDPI1
                     TIJ(2) = DV1DF1*DF1PI2 + DV1DR*DRDPI2 + DV2DF2*DF2PI2 + DV2DR*DRDPI2
                     TIJ(3) = DV1DF1*DF1PI3 + DV1DR*DRDPI3 + DV2DF2*DF2PI3 + DV2DR*DRDPI3

                     DRDPJ1 = DOT_PRODUCT(NR,MATMUL(DRMJ1,-PST(J,:)))
                     DRDPJ2 = DOT_PRODUCT(NR,MATMUL(DRMJ2,-PST(J,:)))
                     DRDPJ3 = DOT_PRODUCT(NR,MATMUL(DRMJ3,-PST(J,:)))

                     TJI(1) = DV1DF1*DF1PJ1 + DV1DR*DRDPJ1 + DV2DF2*DF2PJ1 + DV2DR*DRDPJ1
                     TJI(2) = DV1DF1*DF1PJ2 + DV1DR*DRDPJ2 + DV2DF2*DF2PJ2 + DV2DR*DRDPJ2
                     TJI(3) = DV1DF1*DF1PJ3 + DV1DR*DRDPJ3 + DV2DF2*DF2PJ3 + DV2DR*DRDPJ3

                     G(J3-2:J3) = G(J3-2:J3) - FIJ
                     G(J4-2:J4) = G(J4-2:J4) + FIJ
                     G(J5-2:J5) = G(J5-2:J5) + TIJ
                     G(J6-2:J6) = G(J6-2:J6) + TJI

                  ENDDO

               ENDDO

!     END INNER LOOP OVER PARTICLES

            ENDDO

!     END OUTER LOOP OVER PARTICLES

         ENDDO

         ENERGY = 4.D0*PYEPSNOT*ENERGY
         G      = 24.D0*PYEPSNOT*G

      ELSE

         ENERGY = 0.D0

         DO J1 = 1, REALNATOMS - 1

            J3      = 3*J1
            J5      = OFFSET + J3
            RI      = X(J3-2:J3)
            P       = X(J5-2:J5)

!     ROTATION MATRIX

            THETA  = DSQRT(DOT_PRODUCT(P,P))

            IF (THETA == 0.D0) THEN

               RMI = I3

            ELSE
            
               THETA2 = THETA * THETA
               CT      = COS(THETA)
               ST      = SIN(THETA)
               E(:,:)  = 0.D0
               E(1,2)  = -P(3)
               E(1,3)  =  P(2)
               E(2,3)  = -P(1)
               E(2,1)  = -E(1,2)
               E(3,1)  = -E(1,3)
               E(3,2)  = -E(2,3)
               E       = E/THETA

               RMI     = I3 + (1.D0-CT)*MATMUL(E,E) + E*ST

            ENDIF

!     BEGIN INNER LOOP OVER PARTICLES

            DO J2 = J1 + 1, REALNATOMS

               J4     = 3*J2
               J6     = OFFSET + J4
               RJ     = X(J4-2:J4)
               P      = X(J6-2:J6)

!     ROTATION MATRIX

               THETA  = DSQRT(DOT_PRODUCT(P,P))
              
               IF (THETA == 0.D0) THEN

                  RMJ = I3

               ELSE

                  THETA2 = THETA * THETA
                  CT      = COS(THETA)
                  ST      = SIN(THETA)
                  E(:,:)  = 0.D0
                  E(1,2)  = -P(3)
                  E(1,3)  =  P(2)
                  E(2,3)  = -P(1)
                  E(2,1)  = -E(1,2)
                  E(3,1)  = -E(1,3)
                  E(3,2)  = -E(2,3)
                  E       = E/THETA

                  RMJ     = I3 + (1.D0-CT)*MATMUL(E,E) + E*ST

               ENDIF
 
               DO I = 1, NPYSITE

!     OBTAIN THE SITE POSITION RELATIVE TO THE COM IN THE SPACE-FIXED FRAME

                  RSTI(:) = MATMUL(RMI,PST(I,:))

                  CALL SITEBF(I, EZRI1, EZRI2)

                  AE1 = MATMUL(RMI,(MATMUL(EZRI1(:,:),(TRANSPOSE(RMI)))))

                  IF (RADIFT) THEN

                     AE2 = MATMUL(RMI,(MATMUL(EZRI2(:,:),(TRANSPOSE(RMI)))))

                  ENDIF

                  DO J = 1, NPYSITE

                     RSTJ(:) = MATMUL(RMJ,PST(J,:))
 
                     CALL SITEBF (J, EZRJ1, EZRJ2)

                     BE1 = MATMUL(RMJ,(MATMUL(EZRJ1(:,:),(TRANSPOSE(RMJ)))))

                     IF (RADIFT) THEN

                        BE2 = MATMUL(RMJ,(MATMUL(EZRJ2(:,:),(TRANSPOSE(RMJ)))))

                     ENDIF

!     CALCULATE SEPARATION

                     RIJ    = RI - RJ + RSTI - RSTJ
                     RIJSQ  = DOT_PRODUCT(RIJ,RIJ)
                     ABSRIJ = DSQRT(RIJSQ)
                     NR     = RIJ / ABSRIJ

!     CALCULATE ECF

                     CALL BRENTMIN (0.D0, 0.51D0, 1.D0, AE1, BE1, RIJ, LAMDAC, FMIN)

                     FCNT1   = - FMIN
                     SRTFI1  = 1.D0 / DSQRT(FCNT1)

                     RHO1   = PYSIGNOT / (ABSRIJ - ABSRIJ*SRTFI1 + PYSIGNOT)
                     RHO1SQ = RHO1*RHO1
                     RHO16  = RHO1SQ*RHO1SQ*RHO1SQ
                     RHO112 = RHO16 * RHO16

                     IF (RADIFT) THEN

                        CALL BRENTMIN (0.D0, 0.51D0, 1.D0, AE2, BE2, RIJ, LAMDAC, FMIN)

                        FCNT2   = - FMIN
                        SRTFI2  = 1.D0 / DSQRT(FCNT2)

                        RHO2   = PYSIGNOT / (ABSRIJ - ABSRIJ*SRTFI2 + PYSIGNOT)
                        RHO2SQ = RHO2*RHO2
                        RHO26  = RHO2SQ*RHO2SQ*RHO2SQ

                     ELSE

                        RHO2   = RHO1
                        RHO26  = RHO16

                     ENDIF

                     ENERGY = ENERGY + RHO112 - RHO26

                  ENDDO

               ENDDO

!     CALCULATE PY POTENTIAL ENERGY

            ENDDO

         ENDDO

         ENERGY = 4.D0*PYEPSNOT*ENERGY

      ENDIF

      END SUBROUTINE MULTISITEPY 

!     --------------------------------------------------------------------------

      SUBROUTINE DEFMSPY(PST)

      USE COMMONS, ONLY: PYA1, NPYSITE  

      IMPLICIT NONE

      DOUBLE PRECISION :: LENGTH, PST(NPYSITE,3)

!      LENGTH   = 2.D0*PYA1(1)

      LENGTH = 0.D0
      PST(1,1) = - 0.25D0*LENGTH
      PST(1,2) = 0.D0
      PST(1,3) = 0.D0

      PST(2,1) = 0.25D0*LENGTH
      PST(2,2) = 0.D0
      PST(2,3) = 0.D0

      END SUBROUTINE DEFMSPY

!     --------------------------------------------------------------------------

      SUBROUTINE SITEBF (K, EZR1, EZR2)

      USE COMMONS, ONLY: PYA1, PYA2
    
      IMPLICIT NONE
      
      INTEGER          :: K 
      DOUBLE PRECISION :: EZR1(3,3), EZR2(3,3)

      EZR1(:,:) = 0.D0
      EZR2(:,:) = 0.D0

      IF (K == 1) THEN

         EZR1(1,1) = 1.D0/(PYA1(1)*PYA1(1))
         EZR1(2,2) = 1.D0/(PYA1(2)*PYA1(2))
         EZR1(3,3) = 1.D0/(PYA1(3)*PYA1(3))

         EZR2(1,1) = 1.D0/(PYA2(1)*PYA2(1))
         EZR2(2,2) = 1.D0/(PYA2(2)*PYA2(2))
         EZR2(3,3) = 1.D0/(PYA2(3)*PYA2(3))

      ELSE

         EZR1(1,1) = 1.D0/(PYA1(2)*PYA1(2))
         EZR1(2,2) = 1.D0/(PYA1(1)*PYA1(1))
         EZR1(3,3) = 1.D0/(PYA1(3)*PYA1(3))

         EZR2(1,1) = 1.D0/(PYA2(2)*PYA2(2))
         EZR2(2,2) = 1.D0/(PYA2(1)*PYA2(1))
         EZR2(3,3) = 1.D0/(PYA2(3)*PYA2(3))

      ENDIF

      END SUBROUTINE SITEBF

!     --------------------------------------------------------------------------

      SUBROUTINE DEFINEPYMULTISITES

      USE COMMONS, ONLY: NPYSITE, RADIFT, NATOMS, NSAVE 
      USE PYMODULE, ONLY : PST, OST, ELLST1, ELLST2, ELLMAT, SITECOORDS

      IMPLICIT NONE

      INTEGER          :: J1, J2
      CHARACTER(LEN=5) :: LABEL
      CHARACTER(LEN=7) :: DUMMYLABEL1
      CHARACTER(LEN=11):: DUMMYLABEL2
!      LENGTH   = 2.D0*PYA1(1)

      OPEN(UNIT=299,FILE="PYSITES.XYZ",STATUS="OLD")
      READ(299,*) NPYSITE 
      READ(299,*) 

! ALLOCATE ARRAYS
      IF(.NOT.ALLOCATED(PST)) ALLOCATE(PST(NPYSITE,3))
      IF(.NOT.ALLOCATED(OST)) ALLOCATE(OST(NPYSITE,3))
      IF(.NOT.ALLOCATED(ELLST1)) ALLOCATE(ELLST1(NPYSITE,3))
      IF(.NOT.ALLOCATED(ELLST2)) ALLOCATE(ELLST2(NPYSITE,3))
      IF(.NOT.ALLOCATED(ELLMAT)) ALLOCATE(ELLMAT(NPYSITE,3,3))
      IF(.NOT.ALLOCATED(SITECOORDS)) ALLOCATE(SITECOORDS(NPYSITE,3)) ! ARRAY HOLDING ALL SITE COORDINATES AND ORIENTATIONS
      DO J1=1,NPYSITE
        READ(299,*) LABEL, PST(J1,1), PST(J1,2), PST(J1,3),DUMMYLABEL1,ELLST1(J1,1),ELLST1(J1,2),ELLST1(J1,3),&
                        & ELLST2(J1,1),ELLST2(J1,2),ELLST2(J1,3),DUMMYLABEL2,OST(J1,1),OST(J1,2),OST(J1,3)
      END DO
        ELLST1(:,:)=ELLST1(:,:)/2       ! REPULSIVE SEMIAXES
        ELLST2(:,:)=ELLST2(:,:)/2       ! ATTRACTIVE SEMIAXES

            RADIFT = .TRUE.  ! WILL HAVE TO GET BACK TO THIS ONCE!!!


      END SUBROUTINE DEFINEPYMULTISITES

!     --------------------------------------------------------------------------

      SUBROUTINE PYSITEORIENTATIONS (K, EZR1R, EZR2R)

      USE COMMONS, ONLY: NPYSITE
      USE PYMODULE, ONLY : PST, OST, ELLST1, ELLST2
   
      IMPLICIT NONE
      
      INTEGER          :: K 
      DOUBLE PRECISION :: EZR1(3,3), EZR2(3,3), EZR1R(3,3), EZR2R(3,3)
      DOUBLE PRECISION :: RMAT(3,3),DRMAT1(3,3),DRMAT2(3,3),DRMAT3(3,3)

      EZR1(:,:) = 0.D0
      EZR2(:,:) = 0.D0

! SF344> RMAT IS THE ROTATION MATRIX OF ONE SITE INTO THE BODY-FIXED FRAME, 
!       WE HAVE ONE RMAT FOR EACH ELLIPSOID IN THE RIGID BODY
!    

            CALL RMDRVT(OST(K,:), RMAT, DRMAT1, DRMAT2, DRMAT3, .FALSE.)


! ORIGINAL MATRIX FOR ANGLE-AXIS PART

         EZR1(1,1) = 1.D0/(ELLST1(K,1)*ELLST1(K,1))
         EZR1(2,2) = 1.D0/(ELLST1(K,2)*ELLST1(K,2))
         EZR1(3,3) = 1.D0/(ELLST1(K,3)*ELLST1(K,3))

         EZR2(1,1) = 1.D0/(ELLST2(K,1)*ELLST2(K,1))
         EZR2(2,2) = 1.D0/(ELLST2(K,2)*ELLST2(K,2))
         EZR2(3,3) = 1.D0/(ELLST2(K,3)*ELLST2(K,3))

! NOW ROTATE

        EZR1R = MATMUL(RMAT,(MATMUL(EZR1(:,:),(TRANSPOSE(RMAT)))))
        EZR2R = MATMUL(RMAT,(MATMUL(EZR2(:,:),(TRANSPOSE(RMAT)))))

!        WRITE(*,*) 'BODY ', K
!        WRITE(*,*) 'ORIGINAL MATRIX: '
!        WRITE(*,*) EZR1(:,:)

!        WRITE(*,*) 'ROTATED MATRIX: '
!        WRITE(*,*) EZR1R(:,:)

      END SUBROUTINE PYSITEORIENTATIONS








!     --------------------------------------------------------------------------

      SUBROUTINE MULTISITEPY2 (X, G, ENERGY, GTEST)

      USE COMMONS, ONLY: NATOMS, NPYSITE, PYSIGNOT, PYEPSNOT, RADIFT, MYUNIT
      USE PYMODULE, ONLY : PST, OST
 
      IMPLICIT NONE

      INTEGER          :: I, J, J1, J2, J3, J4, J5, J6, K1, K2, REALNATOMS, OFFSET
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS)
      DOUBLE PRECISION :: EZRI1(3,3), EZRI2(3,3), EZRJ1(3,3), EZRJ2(3,3)
      DOUBLE PRECISION :: RI(3), RJ(3), RIJ(3), NR(3), RIJSQ, ABSRIJ, P(3), THETA, THETA2, CT, ST
      DOUBLE PRECISION :: I3(3,3), AE1(3,3), BE1(3,3), AE2(3,3), BE2(3,3), APB(3,3), APBINV(3,3)
      DOUBLE PRECISION :: RSTI(3), RSTJ(3), FCNT1, FCNT2, SRTFI1, SRTFI2, FMIN, LAMDAC, ENERGY
      DOUBLE PRECISION :: RHO1, RHO1SQ, RHO16, RHO112, RHO2, RHO2SQ, RHO26
      DOUBLE PRECISION :: FCTR1, FCTR2, DV1DF1, DV2DF2, DV1DR, DV2DR
      DOUBLE PRECISION :: DRDPI1, DRDPI2, DRDPI3, DRDPJ1, DRDPJ2, DRDPJ3 
      DOUBLE PRECISION :: DF1PI1, DF1PI2, DF1PI3, DF2PI1, DF2PI2, DF2PI3
      DOUBLE PRECISION :: DF1PJ1, DF1PJ2, DF1PJ3, DF2PJ1, DF2PJ2, DF2PJ3 
      DOUBLE PRECISION :: RMI(3,3), RMJ(3,3), E(3,3)
      DOUBLE PRECISION :: DRMI1(3,3), DRMI2(3,3), DRMI3(3,3), DRMJ1(3,3), DRMJ2(3,3), DRMJ3(3,3)
      DOUBLE PRECISION :: DF1DR(3), DF2DR(3), DG1DR(3), DG2DR(3)
      DOUBLE PRECISION :: ARIBRJ(3), XC(3), XCMRI(3), XCMRJ(3), FIJ(3), TIJ(3), TJI(3)
      DOUBLE PRECISION :: D1ABEZ(3,3), D2ABEZ(3,3), D3ABEZ(3,3), D1ABE(3,3), D2ABE(3,3), D3ABE(3,3) 
      LOGICAL          :: GTEST

      I3(:,:)      = 0.D0
!       NPYSITE = 100
!     FROM INPUT PARAMETERS

      REALNATOMS = NATOMS/2
      OFFSET     = 3*REALNATOMS

      DO K1 = 1, 3

         I3(K1,K1)      = 1.D0

      ENDDO

!      CALL DEFMSPY(PST)

      IF (GTEST) THEN

         ENERGY = 0.D0
         G(:)   = 0.D0

         DO J1 = 1, REALNATOMS - 1

            J3      = 3*J1
            J5      = OFFSET + J3
            RI      = X(J3-2:J3)
            P       = X(J5-2:J5)

!     ROTATION MATRIX

            CALL RMDRVT(P, RMI, DRMI1, DRMI2, DRMI3, GTEST)

!     BEGIN INNER LOOP OVER PARTICLES

            DO J2 = J1 + 1, REALNATOMS

               J4     = 3*J2
               J6     = OFFSET + J4
               RJ     = X(J4-2:J4) 
               P      = X(J6-2:J6)

!     ROTATION MATRIX

               CALL RMDRVT(P, RMJ, DRMJ1, DRMJ2, DRMJ3, GTEST)

               DO I = 1, NPYSITE

!     OBTAIN THE SITE POSITION RELATIVE TO THE COM IN THE SPACE-FIXED FRAME

                  RSTI(:) = MATMUL(RMI,PST(I,:))
                  
                  CALL PYSITEORIENTATIONS(I, EZRI1, EZRI2)

                  AE1 = MATMUL(RMI,(MATMUL(EZRI1(:,:),(TRANSPOSE(RMI)))))

                  IF (RADIFT) THEN

                     AE2 = MATMUL(RMI,(MATMUL(EZRI2(:,:),(TRANSPOSE(RMI)))))

                  ENDIF
               
                  DO J = 1, NPYSITE

                     RSTJ(:) = MATMUL(RMJ,PST(J,:))

                     CALL PYSITEORIENTATIONS(J, EZRJ1, EZRJ2)

                     BE1 = MATMUL(RMJ,(MATMUL(EZRJ1(:,:),(TRANSPOSE(RMJ)))))

                     IF (RADIFT) THEN
   
                        BE2 = MATMUL(RMJ,(MATMUL(EZRJ2(:,:),(TRANSPOSE(RMJ)))))

                     ENDIF

!     CALCULATE SEPARATION

                     RIJ    = RI - RJ + RSTI - RSTJ
                     RIJSQ  = DOT_PRODUCT(RIJ,RIJ)
                     ABSRIJ = DSQRT(RIJSQ)
                     NR     = RIJ / ABSRIJ

!     CALCULATE ECF
                     CALL BRENTMIN (0.D0, 0.51D0, 1.D0, AE1, BE1, RIJ, LAMDAC, FMIN)

                     FCNT1   = - FMIN
                     SRTFI1  = 1.D0 / DSQRT(FCNT1)
                     APB     = LAMDAC * AE1 + (1.D0 - LAMDAC) * BE1

                     CALL MTRXIN (APB, APBINV)

                     ARIBRJ =  LAMDAC * MATMUL(AE1,(RI+RSTI)) + (1.D0 - LAMDAC) * MATMUL(BE1,(RJ+RSTJ))
                     XC     =  MATMUL(APBINV, ARIBRJ)
                     XCMRI  = XC - RI - RSTI
                     XCMRJ  = XC - RJ - RSTJ
                     DF1DR  = - 2.D0 * LAMDAC * MATMUL(AE1,XCMRI)

                     D1ABEZ = MATMUL(DRMI1,EZRI1(:,:))
                     D1ABE  = MATMUL(D1ABEZ,TRANSPOSE(RMI)) + MATMUL(RMI,(TRANSPOSE(D1ABEZ)))

                     D2ABEZ = MATMUL(DRMI2,EZRI1(:,:))
                     D2ABE  = MATMUL(D2ABEZ,(TRANSPOSE(RMI))) + MATMUL(RMI,(TRANSPOSE(D2ABEZ)))

                     D3ABEZ = MATMUL(DRMI3,EZRI1(:,:))
                     D3ABE  = MATMUL(D3ABEZ,(TRANSPOSE(RMI))) + MATMUL(RMI,(TRANSPOSE(D3ABEZ)))

                     DF1PI1 = LAMDAC*(DOT_PRODUCT(XCMRI,MATMUL(D1ABE,XCMRI)) &
                            - 2.D0*DOT_PRODUCT(XCMRI,MATMUL(AE1,MATMUL(DRMI1,PST(I,:))))) 
                     DF1PI2 = LAMDAC*(DOT_PRODUCT(XCMRI,MATMUL(D2ABE,XCMRI)) &
                            - 2.D0*DOT_PRODUCT(XCMRI,MATMUL(AE1,MATMUL(DRMI2,PST(I,:)))))
                     DF1PI3 = LAMDAC*(DOT_PRODUCT(XCMRI,MATMUL(D3ABE,XCMRI)) &
                            - 2.D0*DOT_PRODUCT(XCMRI,MATMUL(AE1,MATMUL(DRMI3,PST(I,:)))))

                     D1ABEZ = MATMUL(DRMJ1,EZRJ1(:,:))
                     D1ABE  = MATMUL(D1ABEZ,(TRANSPOSE(RMJ))) + MATMUL(RMJ,(TRANSPOSE(D1ABEZ)))

                     D2ABEZ = MATMUL(DRMJ2,EZRJ1(:,:))
                     D2ABE  = MATMUL(D2ABEZ,(TRANSPOSE(RMJ))) + MATMUL(RMJ,(TRANSPOSE(D2ABEZ)))

                     D3ABEZ = MATMUL(DRMJ3,EZRJ1(:,:))
                     D3ABE  = MATMUL(D3ABEZ,(TRANSPOSE(RMJ))) + MATMUL(RMJ,(TRANSPOSE(D3ABEZ))) 
               
                     DF1PJ1 = (1.D0-LAMDAC)*(DOT_PRODUCT(XCMRJ,MATMUL(D1ABE,XCMRJ)) &
                            - 2.D0*DOT_PRODUCT(XCMRJ,MATMUL(BE1,MATMUL(DRMJ1,PST(J,:)))))
                     DF1PJ2 = (1.D0-LAMDAC)*(DOT_PRODUCT(XCMRJ,MATMUL(D2ABE,XCMRJ)) &
                            - 2.D0*DOT_PRODUCT(XCMRJ,MATMUL(BE1,MATMUL(DRMJ2,PST(J,:)))))
                     DF1PJ3 = (1.D0-LAMDAC)*(DOT_PRODUCT(XCMRJ,MATMUL(D3ABE,XCMRJ)) &
                            - 2.D0*DOT_PRODUCT(XCMRJ,MATMUL(BE1,MATMUL(DRMJ3,PST(J,:))))) 

                     RHO1   = PYSIGNOT / (ABSRIJ - ABSRIJ*SRTFI1 + PYSIGNOT)
                     RHO1SQ = RHO1*RHO1
                     RHO16  = RHO1SQ*RHO1SQ*RHO1SQ
                     RHO112 = RHO16 * RHO16

                     FCTR1  = 0.5D0*ABSRIJ*SRTFI1/(FCNT1*PYSIGNOT)
                     DG1DR  = (1.D0-SRTFI1)*NR/PYSIGNOT + FCTR1*DF1DR
                     DV1DF1 = -2.D0*RHO112*RHO1*FCTR1
                     DV1DR  = -2.D0*RHO112*RHO1*(1.D0-SRTFI1)/PYSIGNOT

                     IF (RADIFT) THEN

                        CALL BRENTMIN (0.D0, 0.51D0, 1.D0, AE2, BE2, RIJ, LAMDAC, FMIN)

                        FCNT2   = - FMIN
                        SRTFI2  = 1.D0 / DSQRT(FCNT2)
                        APB     = LAMDAC * AE2 + (1.D0 - LAMDAC) * BE2

                        CALL MTRXIN (APB, APBINV)

                        ARIBRJ = LAMDAC * MATMUL(AE2,(RI+RSTI)) + (1.D0 - LAMDAC) * MATMUL(BE2,(RJ+RSTJ))
                        XC     = MATMUL(APBINV, ARIBRJ)
                        XCMRI  = XC - RI - RSTI
                        XCMRJ  = XC - RJ - RSTJ
                        DF2DR  = - 2.D0 * LAMDAC * MATMUL(AE2,XCMRI)

                        RHO2   = PYSIGNOT / (ABSRIJ - ABSRIJ*SRTFI2 + PYSIGNOT)
                        RHO2SQ = RHO2*RHO2
                        RHO26  = RHO2SQ*RHO2SQ*RHO2SQ
               
                        FCTR2  = 0.5D0*ABSRIJ*SRTFI2/(FCNT2*PYSIGNOT)
                        DG2DR  = (1.D0-SRTFI2)*NR/PYSIGNOT+FCTR2*DF2DR
                        DV2DF2 = RHO26*RHO2*FCTR2
                        DV2DR  = RHO26*RHO2*(1.D0-SRTFI2)/PYSIGNOT

                        D1ABEZ = MATMUL(DRMI1,EZRI2(:,:))
                        D1ABE  = MATMUL(D1ABEZ,(TRANSPOSE(RMI))) + MATMUL(RMI,(TRANSPOSE(D1ABEZ)))

                        D2ABEZ = MATMUL(DRMI2,EZRI2(:,:))
                        D2ABE  = MATMUL(D2ABEZ,(TRANSPOSE(RMI))) + MATMUL(RMI,(TRANSPOSE(D2ABEZ)))

                        D3ABEZ = MATMUL(DRMI3,EZRI2(:,:))
                        D3ABE  = MATMUL(D3ABEZ,(TRANSPOSE(RMI))) + MATMUL(RMI,(TRANSPOSE(D3ABEZ)))

                        DF2PI1 = LAMDAC*(DOT_PRODUCT(XCMRI,MATMUL(D1ABE,XCMRI)) &
                               - 2.D0*DOT_PRODUCT(XCMRI,MATMUL(AE2,MATMUL(DRMI1,PST(I,:)))))
                        DF2PI2 = LAMDAC*(DOT_PRODUCT(XCMRI,MATMUL(D2ABE,XCMRI)) &
                               - 2.D0*DOT_PRODUCT(XCMRI,MATMUL(AE2,MATMUL(DRMI2,PST(I,:)))))
                        DF2PI3 = LAMDAC*(DOT_PRODUCT(XCMRI,MATMUL(D3ABE,XCMRI)) &
                               - 2.D0*DOT_PRODUCT(XCMRI,MATMUL(AE2,MATMUL(DRMI3,PST(I,:)))))
                        
                        D1ABEZ = MATMUL(DRMJ1,EZRJ2(:,:))
                        D1ABE  = MATMUL(D1ABEZ,(TRANSPOSE(RMJ))) + MATMUL(RMJ,(TRANSPOSE(D1ABEZ)))

                        D2ABEZ = MATMUL(DRMJ2,EZRJ2(:,:))
                        D2ABE  = MATMUL(D2ABEZ,(TRANSPOSE(RMJ))) + MATMUL(RMJ,(TRANSPOSE(D2ABEZ)))

                        D3ABEZ = MATMUL(DRMJ3,EZRJ2(:,:))
                        D3ABE  = MATMUL(D3ABEZ,(TRANSPOSE(RMJ))) + MATMUL(RMJ,(TRANSPOSE(D3ABEZ)))

                        DF2PJ1 = (1.D0-LAMDAC)*(DOT_PRODUCT(XCMRJ,MATMUL(D1ABE,XCMRJ)) &
                               - 2.D0*DOT_PRODUCT(XCMRJ,MATMUL(BE2,MATMUL(DRMJ1,PST(J,:)))))
                        DF2PJ2 = (1.D0-LAMDAC)*(DOT_PRODUCT(XCMRJ,MATMUL(D2ABE,XCMRJ)) &
                               - 2.D0*DOT_PRODUCT(XCMRJ,MATMUL(BE2,MATMUL(DRMJ2,PST(J,:)))))
                        DF2PJ3 = (1.D0-LAMDAC)*(DOT_PRODUCT(XCMRJ,MATMUL(D3ABE,XCMRJ)) &
                               - 2.D0*DOT_PRODUCT(XCMRJ,MATMUL(BE2,MATMUL(DRMJ3,PST(J,:)))))

                     ELSE

                        RHO2   = RHO1
                        RHO26  = RHO16
                        DG2DR  = DG1DR
                        DV2DF2 = RHO26*RHO2*FCTR1
                        DV2DR  = RHO26*RHO2*(1.D0-SRTFI1)/PYSIGNOT
                        DF2PI1 = DF1PI1
                        DF2PI2 = DF1PI2
                        DF2PI3 = DF1PI3
                        DF2PJ1 = DF1PJ1
                        DF2PJ2 = DF1PJ2
                        DF2PJ3 = DF1PJ3

                     ENDIF
             
!     CALCULATE PY POTENTIAL ENERGY

                     ENERGY = ENERGY + RHO112 - RHO26 
                     IF(ENERGY.LT.-1.0D10) THEN
                        WRITE(MYUNIT,*) 'MULTISITEPY> COLD FUSION DETECTED'
                        ENERGY=1.0D20
                        G=1.0D0
                        RETURN 
                     END IF

!     CALCULATE GRADIENT

                     FIJ = 2.D0*RHO112*RHO1*DG1DR - RHO26*RHO2*DG2DR
                    
                     DRDPI1 = DOT_PRODUCT(NR,MATMUL(DRMI1,PST(I,:)))
                     DRDPI2 = DOT_PRODUCT(NR,MATMUL(DRMI2,PST(I,:)))
                     DRDPI3 = DOT_PRODUCT(NR,MATMUL(DRMI3,PST(I,:)))

                     TIJ(1) = DV1DF1*DF1PI1 + DV1DR*DRDPI1 + DV2DF2*DF2PI1 + DV2DR*DRDPI1
                     TIJ(2) = DV1DF1*DF1PI2 + DV1DR*DRDPI2 + DV2DF2*DF2PI2 + DV2DR*DRDPI2
                     TIJ(3) = DV1DF1*DF1PI3 + DV1DR*DRDPI3 + DV2DF2*DF2PI3 + DV2DR*DRDPI3

                     DRDPJ1 = DOT_PRODUCT(NR,MATMUL(DRMJ1,-PST(J,:)))
                     DRDPJ2 = DOT_PRODUCT(NR,MATMUL(DRMJ2,-PST(J,:)))
                     DRDPJ3 = DOT_PRODUCT(NR,MATMUL(DRMJ3,-PST(J,:)))

                     TJI(1) = DV1DF1*DF1PJ1 + DV1DR*DRDPJ1 + DV2DF2*DF2PJ1 + DV2DR*DRDPJ1
                     TJI(2) = DV1DF1*DF1PJ2 + DV1DR*DRDPJ2 + DV2DF2*DF2PJ2 + DV2DR*DRDPJ2
                     TJI(3) = DV1DF1*DF1PJ3 + DV1DR*DRDPJ3 + DV2DF2*DF2PJ3 + DV2DR*DRDPJ3

                     G(J3-2:J3) = G(J3-2:J3) - FIJ
                     G(J4-2:J4) = G(J4-2:J4) + FIJ
                     G(J5-2:J5) = G(J5-2:J5) + TIJ
                     G(J6-2:J6) = G(J6-2:J6) + TJI

                  ENDDO

               ENDDO

!     END INNER LOOP OVER PARTICLES

            ENDDO

!     END OUTER LOOP OVER PARTICLES

         ENDDO

         ENERGY = 4.D0*PYEPSNOT*ENERGY
         G      = 24.D0*PYEPSNOT*G

      ELSE

         ENERGY = 0.D0

         DO J1 = 1, REALNATOMS - 1

            J3      = 3*J1
            J5      = OFFSET + J3
            RI      = X(J3-2:J3)
            P       = X(J5-2:J5)

!     ROTATION MATRIX

            THETA  = DSQRT(DOT_PRODUCT(P,P))

            IF (THETA == 0.D0) THEN

               RMI = I3

            ELSE
            
               THETA2 = THETA * THETA
               CT      = COS(THETA)
               ST      = SIN(THETA)
               E(:,:)  = 0.D0
               E(1,2)  = -P(3)
               E(1,3)  =  P(2)
               E(2,3)  = -P(1)
               E(2,1)  = -E(1,2)
               E(3,1)  = -E(1,3)
               E(3,2)  = -E(2,3)
               E       = E/THETA

               RMI     = I3 + (1.D0-CT)*MATMUL(E,E) + E*ST

            ENDIF

!     BEGIN INNER LOOP OVER PARTICLES

            DO J2 = J1 + 1, REALNATOMS

               J4     = 3*J2
               J6     = OFFSET + J4
               RJ     = X(J4-2:J4)
               P      = X(J6-2:J6)

!     ROTATION MATRIX

               THETA  = DSQRT(DOT_PRODUCT(P,P))
              
               IF (THETA == 0.D0) THEN

                  RMJ = I3

               ELSE

                  THETA2 = THETA * THETA
                  CT      = COS(THETA)
                  ST      = SIN(THETA)
                  E(:,:)  = 0.D0
                  E(1,2)  = -P(3)
                  E(1,3)  =  P(2)
                  E(2,3)  = -P(1)
                  E(2,1)  = -E(1,2)
                  E(3,1)  = -E(1,3)
                  E(3,2)  = -E(2,3)
                  E       = E/THETA

                  RMJ     = I3 + (1.D0-CT)*MATMUL(E,E) + E*ST

               ENDIF
 
               DO I = 1, NPYSITE

!     OBTAIN THE SITE POSITION RELATIVE TO THE COM IN THE SPACE-FIXED FRAME

                  RSTI(:) = MATMUL(RMI,PST(I,:))

                  CALL SITEBF(I, EZRI1, EZRI2)

                  AE1 = MATMUL(RMI,(MATMUL(EZRI1(:,:),(TRANSPOSE(RMI)))))

                  IF (RADIFT) THEN

                     AE2 = MATMUL(RMI,(MATMUL(EZRI2(:,:),(TRANSPOSE(RMI)))))

                  ENDIF

                  DO J = 1, NPYSITE

                     RSTJ(:) = MATMUL(RMJ,PST(J,:))
 
                     CALL SITEBF (J, EZRJ1, EZRJ2)

                     BE1 = MATMUL(RMJ,(MATMUL(EZRJ1(:,:),(TRANSPOSE(RMJ)))))

                     IF (RADIFT) THEN

                        BE2 = MATMUL(RMJ,(MATMUL(EZRJ2(:,:),(TRANSPOSE(RMJ)))))

                     ENDIF

!     CALCULATE SEPARATION

                     RIJ    = RI - RJ + RSTI - RSTJ
                     RIJSQ  = DOT_PRODUCT(RIJ,RIJ)
                     ABSRIJ = DSQRT(RIJSQ)
                     NR     = RIJ / ABSRIJ

!     CALCULATE ECF

                     CALL BRENTMIN (0.D0, 0.51D0, 1.D0, AE1, BE1, RIJ, LAMDAC, FMIN)

                     FCNT1   = - FMIN
                     SRTFI1  = 1.D0 / DSQRT(FCNT1)

                     RHO1   = PYSIGNOT / (ABSRIJ - ABSRIJ*SRTFI1 + PYSIGNOT)
                     RHO1SQ = RHO1*RHO1
                     RHO16  = RHO1SQ*RHO1SQ*RHO1SQ
                     RHO112 = RHO16 * RHO16

                     IF (RADIFT) THEN

                        CALL BRENTMIN (0.D0, 0.51D0, 1.D0, AE2, BE2, RIJ, LAMDAC, FMIN)

                        FCNT2   = - FMIN
                        SRTFI2  = 1.D0 / DSQRT(FCNT2)

                        RHO2   = PYSIGNOT / (ABSRIJ - ABSRIJ*SRTFI2 + PYSIGNOT)
                        RHO2SQ = RHO2*RHO2
                        RHO26  = RHO2SQ*RHO2SQ*RHO2SQ

                     ELSE

                        RHO2   = RHO1
                        RHO26  = RHO16

                     ENDIF

                     ENERGY = ENERGY + RHO112 - RHO26

                  ENDDO

               ENDDO

!     CALCULATE PY POTENTIAL ENERGY

            ENDDO

         ENDDO

         ENERGY = 4.D0*PYEPSNOT*ENERGY

      ENDIF

      END SUBROUTINE MULTISITEPY2 

!     --------------------------------------------------------------------------

      SUBROUTINE AATOSITES(X,P,XS)
      USE COMMONS, ONLY : NPYSITE
      USE PYMODULE, ONLY : PST,OST,ELLST1,ELLMAT
      IMPLICIT NONE
      
      DOUBLE PRECISION :: X(3),P(3),XS(NPYSITE,3),P1(3),I3(3,3),PHI,THETA,PSI
      DOUBLE PRECISION :: RMAT(3,3), RMAT1(3,3), DRM1(3,3), DRM2(3,3), DRM3(3,3), EZR1(3,3), EZR1R(3,3)
      INTEGER          :: J1
        
        CALL RMDRVT(P, RMAT, DRM1, DRM2, DRM3, .FALSE.)
      I3(:,:)=0.0D0
      DO J1=1,3
        I3(J1,J1)=1.0D0
      END DO
      DO J1=1,NPYSITE
        CALL RMDRVT(OST(J1,:), RMAT1, DRM1, DRM2, DRM3, .FALSE.)

        XS(J1,:)=MATMUL(RMAT,PST(J1,:))
        XS(J1,1)=XS(J1,1) + X(1)
        XS(J1,2)=XS(J1,2) + X(2)
        XS(J1,3)=XS(J1,3) + X(3)

!        AE1 = MATMUL(RMAT,(MATMUL(EZRI1(:,:),(TRANSPOSE(RMAT)))))
! ORIGINAL MATRIX FOR ANGLE-AXIS PART

         EZR1(1,1) = 1.D0/(ELLST1(J1,1)*ELLST1(J1,1))
         EZR1(2,2) = 1.D0/(ELLST1(J1,2)*ELLST1(J1,2))
         EZR1(3,3) = 1.D0/(ELLST1(J1,3)*ELLST1(J1,3))

! NOW ROTATE

        EZR1R = MATMUL(RMAT1,(MATMUL(I3,(TRANSPOSE(RMAT1)))))

! NOW ROTATE AGAIN

!        ELLMAT(J1,:,:) = MATMUL(RMAT,(MATMUL(EZR1R,(TRANSPOSE(RMAT)))))
        ELLMAT(J1,:,:) = MATMUL(RMAT,RMAT1)
!        WRITE(*,*) 'BODY ', J1
!        WRITE(*,*) 'ORIGINAL MATRIX: '
!        WRITE(*,*) EZR1(:,:)

!        WRITE(*,*) 'ROTATION MATRIX: '
!        WRITE(*,*) ELLMAT(J1,:,:)


!   XS: CARTESIAN COORDINATES OF EACH SITE
!   P: ANGLE-AXIS COORDINATES OF THE WHOLE BODY
!        PS(J1,:)=MATMUL(RMAT,OST(J1,:))
      END DO

      END SUBROUTINE AATOSITES

!     --------------------------------------------------------------------------

      SUBROUTINE TAKESTEPMULTISITEPY (NP)

      USE COMMONS
      USE PYMODULE, ONLY : PST, OST
 
      IMPLICIT NONE

      INTEGER          :: I, J, J1, J2, J3, J4, J5, J6, K1, K2, REALNATOMS, OFFSET, NP
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS)
      DOUBLE PRECISION :: EZRI1(3,3), EZRI2(3,3), EZRJ1(3,3), EZRJ2(3,3)
      DOUBLE PRECISION :: RI(3), RJ(3), RIJ(3), NR(3), RIJSQ, ABSRIJ, P(3), THETA, THETA2, CT, ST
      DOUBLE PRECISION :: I3(3,3), AE1(3,3), BE1(3,3), AE2(3,3), BE2(3,3), APB(3,3), APBINV(3,3)
      DOUBLE PRECISION :: RSTI(3), RSTJ(3), FCNT1, FCNT2, SRTFI1, SRTFI2, FMIN, LAMDAC, ENERGY
      DOUBLE PRECISION :: RHO1, RHO1SQ, RHO16, RHO112, RHO2, RHO2SQ, RHO26
      DOUBLE PRECISION :: FCTR1, FCTR2, DV1DF1, DV2DF2, DV1DR, DV2DR
      DOUBLE PRECISION :: DRDPI1, DRDPI2, DRDPI3, DRDPJ1, DRDPJ2, DRDPJ3 
      DOUBLE PRECISION :: DF1PI1, DF1PI2, DF1PI3, DF2PI1, DF2PI2, DF2PI3
      DOUBLE PRECISION :: DF1PJ1, DF1PJ2, DF1PJ3, DF2PJ1, DF2PJ2, DF2PJ3 
      DOUBLE PRECISION :: RMI(3,3), RMJ(3,3), E(3,3)
      DOUBLE PRECISION :: DRMI1(3,3), DRMI2(3,3), DRMI3(3,3), DRMJ1(3,3), DRMJ2(3,3), DRMJ3(3,3)
      DOUBLE PRECISION :: DF1DR(3), DF2DR(3), DG1DR(3), DG2DR(3)
      DOUBLE PRECISION :: ARIBRJ(3), XC(3), XCMRI(3), XCMRJ(3), FIJ(3), TIJ(3), TJI(3)
      DOUBLE PRECISION :: D1ABEZ(3,3), D2ABEZ(3,3), D3ABEZ(3,3), D1ABE(3,3), D2ABE(3,3), D3ABE(3,3) 
      DOUBLE PRECISION :: RANDOM,DPRAND,ECFVAL,LOCALSTEP,DUMMY2 
      LOGICAL          :: OVERLAPT

      I3(:,:)      = 0.D0

      REALNATOMS = NATOMS/2
      OFFSET     = 3*REALNATOMS

      DO K1 = 1, 3

         I3(K1,K1)      = 1.D0

      ENDDO


         DO J1 = 1, REALNATOMS - 1


OVERLAPT=.TRUE.

95 DO WHILE(OVERLAPT)
            LOCALSTEP = 0.0D0
            J3      = 3*J1
            J5      = OFFSET + J3

         DUMMY2 = COORDS(J3-2,NP)**2 + COORDS(J3-1,NP)**2 + COORDS(J3,NP)**2
      
         IF (DUMMY2 .GT. RADIUS) THEN
! BRING BACK THE MOLECULE WITHIN THE RADIUS
                COORDS(J3-2,NP)=COORDS(J3-2,NP)-SQRT(RADIUS)*NINT(COORDS(J3-2,NP)/SQRT(RADIUS))
                COORDS(J3-1,NP)=COORDS(J3-1,NP)-SQRT(RADIUS)*NINT(COORDS(J3-1,NP)/SQRT(RADIUS))
                COORDS(J3,NP)=COORDS(J3,NP)-SQRT(RADIUS)*NINT(COORDS(J3,NP)/SQRT(RADIUS))
    WRITE(MYUNIT,'(A,2F20.10)') 'INITIAL COORDINATE OUTSIDE CONTAINER -- BRINGING MOLECULE BACK WITHIN THE CONTAINER RADIUS'
        END IF
!     BEGIN INNER LOOP OVER PARTICLES
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


            IF(FROZEN(J1)) LOCALSTEP = 0.0D0

            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
            COORDS(J5-2,NP) = COORDS(J5-2,NP) + LOCALSTEP*RANDOM
            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
            COORDS(J5-1,NP) = COORDS(J5-1,NP) + LOCALSTEP*RANDOM
            RANDOM          = (DPRAND() - 0.5D0)*2.0D0
            COORDS(J5,NP)   = COORDS(J5,NP) + LOCALSTEP*RANDOM
!   
            OVERLAPT = .FALSE.

            X(J3-2:J3) = COORDS(J3-2:J3,NP)
            X(J5-2:J5) = COORDS(J5-2:J5,NP)
            RI      = X(J3-2:J3)
            P       = X(J5-2:J5)

!     ROTATION MATRIX

            CALL RMDRVT(P, RMI, DRMI1, DRMI2, DRMI3, .FALSE.)


            DO J2 = J1 + 1, REALNATOMS

               J4     = 3*J2
               J6     = OFFSET + J4
               X(J4-2:J4) = COORDS(J4-2:J4,NP)
               X(J6-2:J6) = COORDS(J6-2:J6,NP)
               RJ     = X(J4-2:J4) 
               P      = X(J6-2:J6)

!     ROTATION MATRIX

               CALL RMDRVT(P, RMJ, DRMJ1, DRMJ2, DRMJ3, .FALSE.)

               DO I = 1, NPYSITE

!     OBTAIN THE SITE POSITION RELATIVE TO THE COM IN THE SPACE-FIXED FRAME

                  RSTI(:) = MATMUL(RMI,PST(I,:))
                  
                  CALL PYSITEORIENTATIONS(I, EZRI1, EZRI2)

                  AE1 = MATMUL(RMI,(MATMUL(EZRI1(:,:),(TRANSPOSE(RMI)))))

                  IF (RADIFT) THEN

                     AE2 = MATMUL(RMI,(MATMUL(EZRI2(:,:),(TRANSPOSE(RMI)))))

                  ENDIF
               
                  DO J = 1, NPYSITE

                     RSTJ(:) = MATMUL(RMJ,PST(J,:))

                     CALL PYSITEORIENTATIONS(J, EZRJ1, EZRJ2)

                     BE1 = MATMUL(RMJ,(MATMUL(EZRJ1(:,:),(TRANSPOSE(RMJ)))))

                     IF (RADIFT) THEN
   
                        BE2 = MATMUL(RMJ,(MATMUL(EZRJ2(:,:),(TRANSPOSE(RMJ)))))

                     ENDIF

!     CALCULATE SEPARATION

                     RIJ    = RI - RJ + RSTI - RSTJ
                     RIJSQ  = DOT_PRODUCT(RIJ,RIJ)
                     ABSRIJ = DSQRT(RIJSQ)
                     NR     = RIJ / ABSRIJ

!     CALCULATE ECF
                     CALL BRENTMIN (0.D0, 0.51D0, 1.D0, AE1, BE1, RIJ, LAMDAC, FMIN)
                          ECFVAL = - FMIN
!                     ALLOW FOR A SLIGHT OVERLAP 
                          IF (ECFVAL < PYOVERLAPTHRESH) THEN
!                                WRITE(*,*) 'ATOMS OVERLAPPING', J1, J2, ECFVAL, ABSRIJ
                             OVERLAPT = .TRUE.
                             GO TO 95
                          ENDIF

            

                  ENDDO  ! END LOOP OVER SITES IN THE SECOND BODY

               ENDDO ! END LOOP OVER SITES IN THE FIRST BODY

!     END INNER LOOP OVER PARTICLES

            ENDDO

         ENDDO ! END WHILE


!     END OUTER LOOP OVER PARTICLES

         ENDDO



      END SUBROUTINE TAKESTEPMULTISITEPY 

!     --------------------------------------------------------------------------


