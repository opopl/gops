
! QO => R{{{
        DO I = 1, N
          J = (I-1)*3
          DO K=1,3
            R(I,K) = QO(J+K)
          ENDDO
        ENDDO
! }}}
! INTER-PARTICLE DISTANCES: DR, LEN_DR; BOND VECTORS: BVR {{{

        DO I = 1, N-1
          DO J = I+1, N
            DR(I,J,1:3) = R(J,1:3) - R(I,1:3)
            DR(J,I,1:3) = -DR(I,J,1:3)
            LEN_DR(I,J) = SQRT(SUM(DR(I,J,1:3)**2))
            LEN_DR(J,I) = LEN_DR(I,J)
          ENDDO
            BVR(I,1:3)=DR(I,I+1,1:3)
        ENDDO
! }}}
! DOT PRODUCTS BETWEEN BOND VECTORS: DPD, LEN_BVR  {{{

      DO I = 1, N-1
        IF ( I .LE. N-3 ) THEN KMAX=3
        IF ( I .EQ. N-2 ) THEN KMAX=2
        IF ( I .LE. N-1 ) THEN KMAX=1
        DO K=1,KMAX
         J=I+K-1
         DPD(I,K) = SUM(BVR(I,1:3)*BVR(J,1:3))
        ENDDO
        LEN_BVR(I)=SQRT(DPD(I,1))
      ENDDO
! }}}
! Squared cross-products between adjacent bond vectors i and i+1: XPD_2 {{{

        DO I = 1, N-2
           XPD_2(I) = DPD(I,1)*DPD(I+1,1)-DPD(I,2)**2 
           XPD(I)=SQRT(XPD_2(I))
           VXPD(I,1)=BVR(I,2)*BVR(I+1,3)-BVR(I,3)*BVR(I+1,2)
           VXPD(I,2)=BVR(I,3)*BVR(I+1,1)-BVR(I,1)*BVR(I+1,3)
           VXPD(I,3)=BVR(I,1)*BVR(I+1,2)-BVR(I,2)*BVR(I+1,1)
           HVXPD(I,1:3)=VXPD(I,1:3)/XPD(I) 
        ENDDO

! }}}
! BOND ANGLES: ANG(I,1), I=2,...,N-1 {{{

        DO I = 1, N-2
            COS_THETA=-DPD(I,2)/(LEN_BVR(I)*LEN_BVR(I+1))
            ANG(I+1,1) = ACOS(COS_THETA)
        ENDDO
! }}}
! TORSIONAL ANGLES: ANG(I,2), I=2,...,N-2 {{{

        DO I = 1, N-3
            COS_PHI = (DPD(I,2)*DPD(I+1,2)-DPD(I,3)*DPD(I+1,1))
            COS_PHI = COS_PHI/(XPD(I)*XPD(I+1))
            IF (ABS(COS_PHI).GT.1.0D0) COS_PHI=COS_PHI/ABS(COS_PHI)
            ANG(I+1,2) = ACOS(COS_PHI)
        ENDDO
! }}}

