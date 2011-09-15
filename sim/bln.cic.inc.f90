
! INTER-PARTICLE DISTANCES: DR, LEN_DR; BOND VECTORS: BVR {{{

        DO I = 1, N-1
          DO J = I+1, N
            DR(I,J,1:3) = R(I,1:3) - R(J,1:3)
            DR(J,I,1:3) = -DR(I,J,1:3)
            LEN_DR(I,J) = SQRT(SUM(DR(I,J,1:3)**2))
            LEN_DR(J,I) = LEN_DR(I,J)
          ENDDO
            BVR(I,1:3)=DR(I+1,I,1:3)
        ENDDO
! }}}
! BOND VECTORS: DPD, B, EB  {{{

      DO I = 1, N-1
        B(I)= SUM(BVR(I,1:3)**2)
        EB(I,1:3)=BVR(I,1:3)/B(I)
        IF (I.LT.N-1) DPD(I)= SUM(BVR(I,1:3)*BVR(I+1,1:3))
      ENDDO
! }}}
! Cross-products between adjacent bond vectors i and i+1 {{{

        DO I = 1, N-2
           XPD_2(I) = B(I)*B(I+1)-DPD(I)**2 
           XPD(I)=SQRT(XPD_2(I))
           VXPD(I,1)=BVR(I,2)*BVR(I+1,3)-BVR(I,3)*BVR(I+1,2)
           VXPD(I,2)=BVR(I,3)*BVR(I+1,1)-BVR(I,1)*BVR(I+1,3)
           VXPD(I,3)=BVR(I,1)*BVR(I+1,2)-BVR(I,2)*BVR(I+1,1)
           HVXPD(I,1:3)=VXPD(I,1:3)/XPD(I) 
           PP(I,1:3)=VXPD(I,1:3)/XPD_2(I)
        ENDDO

! }}}
! BOND ANGLES: ANG(I,1), I=2,...,N-1 {{{

        DO I = 2, N-1
            COS_THETA=-DPD(I-1)/(B(I-1)*B(I))
            ANG(I,1) = ACOS(COS_THETA)
            F(I,1)=2.0D0*RK_THETA*(ANG(I,1)-THETA_0)
        ENDDO
! }}}
! TORSIONAL ANGLES: ANG(I,2), I=2,...,N-2 {{{

        DO I = 2, N-2
            COS_PHI=-SUM(HVXPD(I-1,1:3)*HVXPD(I,1:3))
            IF (ABS(COS_PHI).GT.1.0D0) COS_PHI=COS_PHI/ABS(COS_PHI)
            ANG(I,2) = ACOS(COS_PHI)
            AN=ANG(I,2)
            F(I,2)=-CD(I,1)*SIN(AN)-3.0*CD(I,2)*SIN(3.0*AN)
            FB(I,1)=F(I,2)*B(I)
            FB(I,2)=F(I,2)/B(I)
        ENDDO
! }}}

