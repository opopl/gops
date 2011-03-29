C   GPL LICENSE INFO {{{
C   GMIN: A PROGRAM FOR FINDING GLOBAL MINIMA
C   COPYRIGHT (C) 1999-2006 DAVID J. WALES
C   THIS FILE IS PART OF GMIN.
C
C   GMIN IS FREE SOFTWARE; YOU CAN REDISTRIBUTE IT AND/OR MODIFY
C   IT UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE AS PUBLISHED BY
C   THE FREE SOFTWARE FOUNDATION; EITHER VERSION 2 OF THE LICENSE, OR
C   (AT YOUR OPTION) ANY LATER VERSION.
C
C   GMIN IS DISTRIBUTED IN THE HOPE THAT IT WILL BE USEFUL,
C   BUT WITHOUT ANY WARRANTY; WITHOUT EVEN THE IMPLIED WARRANTY OF
C   MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  SEE THE
C   GNU GENERAL PUBLIC LICENSE FOR MORE DETAILS.
C
C   YOU SHOULD HAVE RECEIVED A COPY OF THE GNU GENERAL PUBLIC LICENSE
C   ALONG WITH THIS PROGRAM; IF NOT, WRITE TO THE FREE SOFTWARE
C   FOUNDATION, INC., 59 TEMPLE PLACE, SUITE 330, BOSTON, MA  02111-1307  USA
C
C }}}

C> CALCULATE THE ENERGY AND GRADIENT
C> FOR A GIVEN CONFIGURATION OF A BLN POLYMER CHAIN.

      SUBROUTINE BLN(QO,GRAD,ENERGY,GRADT)
C{{{
      USE COMMONS
      IMPLICIT NONE
      DOUBLE PRECISION QO(3*NATOMS), GRAD(3*NATOMS)
      LOGICAL GRADT
      INTEGER N
      DOUBLE PRECISION X(NATOMS), Y(NATOMS), Z(NATOMS), XR(NATOMS,NATOMS), YR(NATOMS,NATOMS), 
     &                 ZR(NATOMS,NATOMS), DOT_PROD(NATOMS,3),
     &                 X_PROD(NATOMS), BOND_ANGLE(NATOMS), TOR_ANGLE(NATOMS), RADII(NATOMS,NATOMS), 
     &                 ENERGY, COSTOR(NATOMS), DFAC(NATOMS), SINBOND(NATOMS)

      N=NATOMS
!
! WITHOUT THESE INITIALISATIONS THE NAG COMPILER FILLS IN RANDOM NUMBERS FOR
! UNASSIGNED ELEMENTS WITH OPTIMISATION TURNED ON.
!
      BOND_ANGLE(1:NATOMS)=0.0D0
      TOR_ANGLE(1:NATOMS)=0.0D0
      COSTOR(1:NATOMS)=0.0D0
      DFAC(1:NATOMS)=0.0D0
      SINBOND(1:NATOMS)=0.0D0
      DOT_PROD(1:NATOMS,1:3)=0.0D0
      X_PROD(1:NATOMS)=0.0D0
      RADII(1:NATOMS,1:NATOMS)=0.0D0
      XR(1:NATOMS,1:NATOMS)=0.0D0
      YR(1:NATOMS,1:NATOMS)=0.0D0
      ZR(1:NATOMS,1:NATOMS)=0.0D0

      CALL CALC_INT_COORDSBLN(QO,N,X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD,BOND_ANGLE,TOR_ANGLE,RADII,NATOMS,COSTOR,
     &                        SINBOND,A_BLN,B_BLN,C_BLN,D_BLN,DFAC)
      CALL CALC_ENERGYBLN(QO,ENERGY,N,LJREP_BLN,LJATT_BLN,A_BLN,B_BLN,C_BLN,D_BLN,X,Y,Z,XR,YR,ZR,DOT_PROD,
     &                    X_PROD,BOND_ANGLE,TOR_ANGLE,RADII,NATOMS,RK_R,RK_THETA,COSTOR)
      IF (.NOT.GRADT) RETURN
 
      CALL CALC_GRADIENTBLN(QO,GRAD,N,LJREP_BLN,LJATT_BLN,A_BLN,B_BLN,C_BLN,D_BLN,X,Y,Z,XR,YR,ZR,DOT_PROD,
     &                      X_PROD,BOND_ANGLE,TOR_ANGLE,RADII,NATOMS,RK_R,RK_THETA,COSTOR,DFAC,SINBOND)

      RETURN
      END
C }}}
C
C> CALCULATE THE INTERNAL COORDINATES
C
      SUBROUTINE CALC_INT_COORDSBLN(QO,N,X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD,BOND_ANGLE,TOR_ANGLE,RADII,NATOMS, 
     &                              COSTOR,SINBOND,A_BLN,B_BLN,C_BLN,D_BLN,DFAC)
C {{{
C DECLARATIONS {{{
      IMPLICIT NONE
      INTEGER I, N, J, NATOMS
      DOUBLE PRECISION COS_PHI, COS_THETA, DUMMY, DUMMY2
      DOUBLE PRECISION QO(3*NATOMS)
      DOUBLE PRECISION, PARAMETER :: TWOPI=6.283185307179586477D0
      DOUBLE PRECISION X(NATOMS), Y(NATOMS), Z(NATOMS), XR(NATOMS,NATOMS), DFAC(NATOMS),
     &  YR(NATOMS,NATOMS), ZR(NATOMS,NATOMS), DOT_PROD(NATOMS,3), COSTOR(NATOMS), D_BLN(NATOMS), 
     &  A_BLN(NATOMS), B_BLN(NATOMS), X_PROD(NATOMS), BOND_ANGLE(NATOMS), TOR_ANGLE(NATOMS), 
     &  RADII(NATOMS,NATOMS), SINBOND(NATOMS), C_BLN(NATOMS)
C }}}
      DO I = 1, N
         J = (I-1)*3
         X(I) = QO((I-1)*3+1)
         Y(I) = QO((I-1)*3+2)
         Z(I) = QO((I-1)*3+3)
      ENDDO
C
C INTER-PARTICLE DISTANCES {{{
C
      DO I = 1, N-1
         DO J = I+1, N
C        DO J = 1, N
            XR(I,J) = X(J) - X(I)
            YR(I,J) = Y(J) - Y(I)
            ZR(I,J) = Z(J) - Z(I)
            RADII(I,J) = SQRT(XR(I,J)*XR(I,J) + YR(I,J)*YR(I,J) + ZR(I,J)*ZR(I,J))
            RADII(J,I) = RADII(I,J)
         ENDDO
      ENDDO
C }}}
C
C DOT PRODUCTS BETWEEN BOND VECTORS {{{
C

      DO I = 1, N-3
         DOT_PROD(I,1) = XR(I,I+1)*XR(I,I+1) + YR(I,I+1)*YR(I,I+1) + ZR(I,I+1)*ZR(I,I+1)
         DOT_PROD(I,2) = XR(I,I+1)*XR(I+1,I+2)+YR(I,I+1)*YR(I+1,I+2)+ ZR(I,I+1)*ZR(I+1,I+2)
         DOT_PROD(I,3) = XR(I,I+1)*XR(I+2,I+3)+YR(I,I+1)*YR(I+2,I+3)+ ZR(I,I+1)*ZR(I+2,I+3)
      ENDDO

!     I = N-2
      DOT_PROD(N-2,1) = XR(N-2,N-2+1)*XR(N-2,N-2+1) + YR(N-2,N-2+1)*YR(N-2,N-2+1) + ZR(N-2,N-2+1)*ZR(N-2,N-2+1)
      DOT_PROD(N-2,2) = XR(N-2,N-2+1)*XR(N-2+1,N-2+2)+YR(N-2,N-2+1)*YR(N-2+1,N-2+2)+ ZR(N-2,N-2+1)*ZR(N-2+1,N-2+2)
!     I = N-1
      DOT_PROD(N-1,1) = XR(N-1,N-1+1)*XR(N-1,N-1+1) + YR(N-1,N-1+1)*YR(N-1,N-1+1) + ZR(N-1,N-1+1)*ZR(N-1,N-1+1)
C }}}
C CROSS-PRODUCTS BETWEEN ADJACENT BOND VECTORS {{{

      DO I = 1, N-2
         X_PROD(I) = DOT_PROD(I,1)*DOT_PROD(I+1,1) - DOT_PROD(I,2)*DOT_PROD(I,2)   
      ENDDO
C }}}
C BOND ANGLES {{{

      DO I = 1, N-2
         COS_THETA=-DOT_PROD(I,2)/(SQRT(DOT_PROD(I,1)*DOT_PROD(I+1,1)))
         BOND_ANGLE(I+1) = DACOS(COS_THETA)
         SINBOND(I+1)=SIN(BOND_ANGLE(I+1))*SQRT(DOT_PROD(I,1)*DOT_PROD(I+1,1))
      ENDDO
C }}}
C TORSIONAL ANGLES {{{

      DO I = 1, N-3
         COS_PHI = (DOT_PROD(I,2)*DOT_PROD(I+1,2) - DOT_PROD(I,3)*DOT_PROD(I+1,1))/SQRT(X_PROD(I)*X_PROD(I+1))
         IF (ABS(COS_PHI).GT.1.0D0) COS_PHI=SIGN(1.0D0,COS_PHI)
         TOR_ANGLE(I+1) = DACOS(COS_PHI)
C }}}
C
C TOR_ANGLE IS RETURNED IN THE RANGE 0 TO PI.
C DUMMY SHOULD TAKE THE OPPOSITE SIGN FROM THE DIHEDRAL ANGLE. NEGATIVE
C VALUES OF THE DIHEDRAL SHOULD BE SUBTRACTED FROM 2*PI.
C THIS IS ONLY NECESSARY WHEN THE POTENTIAL CONTAINS COS(PHI+PI/4) TERMS BECAUSE
C THE GRADIENT IS DISCONTINUOUS WHEN PHI GOES THROUGH PI FOR SUCH TERMS IF PHI
C IS RESTRICTED TO 0 < PHI < PI.
C {{{
C        DUMMY=XR(I+2,I+3)*(YR(I+1,I)*ZR(I+1,I+2)-YR(I+1,I+2)*ZR(I+1,I))+
C    &         YR(I+2,I+3)*(XR(I+1,I+2)*ZR(I+1,I)-XR(I+1,I)*ZR(I+1,I+2))+
C    &         ZR(I+2,I+3)*(XR(I+1,I)*YR(I+1,I+2)-XR(I+1,I+2)*YR(I+1,I))
         DUMMY=XR(I+2,I+3)*(-YR(I,I+1)*ZR(I+1,I+2)+YR(I+1,I+2)*ZR(I,I+1))+
     &         YR(I+2,I+3)*(-XR(I+1,I+2)*ZR(I,I+1)+XR(I,I+1)*ZR(I+1,I+2))+
     &         ZR(I+2,I+3)*(-XR(I,I+1)*YR(I+1,I+2)+XR(I+1,I+2)*YR(I,I+1))
         IF (DUMMY.GT.0.0D0) TOR_ANGLE(I+1)=TWOPI-TOR_ANGLE(I+1)
         COSTOR(I+1)=COS(TOR_ANGLE(I+1))
C }}}
C  THIS IS AN UGLY HACK TO PREVENT DIVISION BY ZERO. THERE WILL BE A LOSS OF PRECISION
C  IF A DIHEDRAL IS 0, PI, 2*PI, IF D_BLN IS NON-ZERO.
C {{{
         IF (TAN(TOR_ANGLE(I+1)).EQ.0.0D0) THEN
            PRINT '(A,I8,A,G20.10)','WARNING IN BLN, DIHEDRAL ANGLE ',I+1,' IS ',TOR_ANGLE(I+1)
            TOR_ANGLE(I+1)=TOR_ANGLE(I+1)+1.0D-10
            PRINT '(A,G20.10)','WARNING IN BLN, TAN PERTURBED ANGLE=',TOR_ANGLE(I+1)
         ENDIF
         DUMMY2=TAN(TOR_ANGLE(I+1))
         DFAC(I+1)=(A_BLN(I+1)+D_BLN(I+1)*( 1.0D0+1.0D0/DUMMY2 )*0.7071067811865475244D0-B_BLN(I+1)
     &              +C_BLN(I+1)*(12.0*COSTOR(I+1)**2-3.0))/SQRT(X_PROD(I+1)*X_PROD(I))
C }}}
      ENDDO

      RETURN
      END
C }}}
C> CALCULATE THE ENERGY

      SUBROUTINE CALC_ENERGYBLN(QO,ENERGY,N,LJREP_BLN,LJATT_BLN,A_BLN,B_BLN,C_BLN,D_BLN,X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD,
     &  BOND_ANGLE,TOR_ANGLE,RADII,NATOMS,RK_R,RK_THETA,COSTOR)
C {{{
      IMPLICIT NONE
      INTEGER NATOMS, N, I, J
      DOUBLE PRECISION, PARAMETER :: THETA_0 = 1.8326D0, PI4=0.7853981633974483096D0 ! 1.8326 RADIANS IS 105 DEGREES
      DOUBLE PRECISION QO(3*NATOMS), ENERGY
      DOUBLE PRECISION X(NATOMS), Y(NATOMS), Z(NATOMS), XR(NATOMS,NATOMS), 
     &  YR(NATOMS,NATOMS), ZR(NATOMS,NATOMS), DOT_PROD(NATOMS,3), COSTOR(NATOMS),
     &  X_PROD(NATOMS), BOND_ANGLE(NATOMS), TOR_ANGLE(NATOMS), RADII(NATOMS,NATOMS)
      DOUBLE PRECISION RK_R, RK_THETA, E_NBOND, E_BOND, E_BANGLE, E_TANGLE, RAD6
      DOUBLE PRECISION LJREP_BLN(NATOMS,NATOMS),LJATT_BLN(NATOMS,NATOMS),C_BLN(NATOMS),A_BLN(NATOMS),B_BLN(NATOMS),D_BLN(NATOMS)

      E_NBOND=0.0D0
      E_BOND=0.0D0
      E_BANGLE=0.0D0
      E_TANGLE=0.0D0

      DO I = 1, N-2
         DO J = I+2, N
            RAD6 = RADII(I,J)**6
            E_NBOND = E_NBOND + (LJREP_BLN(I,J)/RAD6 + LJATT_BLN(I,J))/RAD6
         ENDDO
      ENDDO
      E_NBOND=E_NBOND*4.0D0

      DO I = 1, N-1
         E_BOND = E_BOND + (RADII(I,I+1)-1.0D0)**2
      ENDDO
      E_BOND=E_BOND*RK_R/2.0D0

      DO I = 2, N-1
         E_BANGLE = E_BANGLE + (BOND_ANGLE(I)-THETA_0)**2
      ENDDO
      E_BANGLE=E_BANGLE*RK_THETA/2.0D0

      DO I = 2, N-2
         E_TANGLE = E_TANGLE + A_BLN(I)*(1.0D0 + COSTOR(I)) + C_BLN(I)*(1.0 + COS(3.0*TOR_ANGLE(I)))
     &                       + B_BLN(I)*(1.0D0 - COSTOR(I)) + D_BLN(I)*(1.0 + COS(TOR_ANGLE(I)+PI4))
      ENDDO

      ENERGY = E_NBOND + E_BOND + E_BANGLE + E_TANGLE

C     WRITE(*,'(A,4F20.10)') 'NBOND,BOND,BANGLE,TANGLE=',E_NBOND,E_BOND,E_BANGLE,E_TANGLE

      RETURN
      END
C }}}
C
C> CALCULATE THE GRADIENTS
C
      SUBROUTINE CALC_GRADIENTBLN(QO,FQ,N,LJREP_BLN,LJATT_BLN,A_BLN,B_BLN,C_BLN,D_BLN,
     &                            X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD,
     &                            BOND_ANGLE,TOR_ANGLE,RADII,NATOMS,RK_R,RK_THETA,COSTOR,DFAC,SINBOND)
C {{{
C DECLARATIONS {{{
      USE COMMONS,ONLY : MYUNIT
      IMPLICIT NONE
      INTEGER NATOMS, N, I, J
      DOUBLE PRECISION, PARAMETER :: THETA_0 = 1.8326D0
      DOUBLE PRECISION QO(3*NATOMS),FQ(3*NATOMS),FX(NATOMS),FY(NATOMS), FZ(NATOMS), DFAC(NATOMS), SINBOND(NATOMS)
      DOUBLE PRECISION FNB_X(NATOMS),FNB_Y(NATOMS),FNB_Z(NATOMS), FB_X(NATOMS),FB_Y(NATOMS)
      DOUBLE PRECISION FB_Z(NATOMS),FBA_X(NATOMS),FBA_Y(NATOMS),FBA_Z(NATOMS), COSTOR(NATOMS)
      DOUBLE PRECISION FTA_X(NATOMS),FTA_Y(NATOMS),FTA_Z(NATOMS), A3, COEF, COEF1, COEF2, COEF3, A4
      DOUBLE PRECISION RK_R, RK_THETA, RAD7, RAD14, DF, FXX, FZZ, FYY, RVAR, DEN, RNUM, DEN1, A1, A2, DEN2

      DOUBLE PRECISION X(NATOMS), Y(NATOMS), Z(NATOMS), XR(NATOMS,NATOMS), 
     &  YR(NATOMS,NATOMS), ZR(NATOMS,NATOMS), DOT_PROD(NATOMS,3), 
     &  X_PROD(NATOMS), BOND_ANGLE(NATOMS), TOR_ANGLE(NATOMS), RADII(NATOMS,NATOMS)

      DOUBLE PRECISION LJREP_BLN(NATOMS,NATOMS), LJATT_BLN(NATOMS,NATOMS),C_BLN(NATOMS),A_BLN(NATOMS),B_BLN(NATOMS),D_BLN(NATOMS)
C }}}
C
C GRADIENTS OF POTENTIAL
C {{{
      DO I = 1,N

         FNB_X(I) = 0.0  
         FNB_Y(I) = 0.0 
         FNB_Z(I) = 0.0 
   
         FB_X(I)  = 0.0 
         FB_Y(I)  = 0.0 
         FB_Z(I)  = 0.0 
   
         FBA_X(I) = 0.0 
         FBA_Y(I) = 0.0 
         FBA_Z(I) = 0.0 
   
         FTA_X(I) = 0.0 
         FTA_Y(I) = 0.0 
         FTA_Z(I) = 0.0 
   
         FX(I)= 0.0 
         FY(I)= 0.0 
         FZ(I)= 0.0 

      ENDDO
C }}}
C ..... NON-BONDED INTERACTION FORCES ..... 
C {{{
      DO I = 1, N-2
         DO J = I+2, N

            RAD7 = RADII(I,J)*RADII(I,J)*RADII(I,J)*RADII(I,J)*RADII(I,J)*RADII(I,J)*RADII(I,J)   
            RAD14 = RAD7*RAD7 

            DF = -24.0*((2.0*LJREP_BLN(I,J)/RAD14) + (LJATT_BLN(I,J)/(RAD7*RADII(I,J))))

            FXX = DF*XR(I,J) 
            FYY = DF*YR(I,J) 
            FZZ = DF*ZR(I,J) 

            FNB_X(I) = FXX + FNB_X(I)
            FNB_Y(I) = FYY + FNB_Y(I)
            FNB_Z(I) = FZZ + FNB_Z(I)

            FNB_X(J) = -FXX + FNB_X(J)
            FNB_Y(J) = -FYY + FNB_Y(J)
            FNB_Z(J) = -FZZ + FNB_Z(J)

         ENDDO
      ENDDO
C }}}
C ... BOND INTERACTION FORCES ... 
C {{{
      DO I = 1, N-1

         RVAR = 1.0D0/RADII(I,I+1) 

         DF = RK_R*(1.0 - RVAR) 
         FXX = DF*XR(I,I+1) 
         FYY = DF*YR(I,I+1) 
         FZZ = DF*ZR(I,I+1) 

         FB_X(I) = FXX + FB_X(I)
         FB_Y(I) = FYY + FB_Y(I)
         FB_Z(I) = FZZ + FB_Z(I)

         FB_X(I+1) = -FXX + FB_X(I+1)
         FB_Y(I+1) = -FYY + FB_Y(I+1)
         FB_Z(I+1) = -FZZ + FB_Z(I+1)

      ENDDO
C }}}
C BOND ANGLE FORCES  PARTICLE 1
C PARTICLES 1,2,N-1, AND N DONE OUTSIDE OF THE LOOP
C {{{
!     I = 1
      DEN = SINBOND(1+1)
      RNUM = RK_THETA*(BOND_ANGLE(1+1) - THETA_0)

      FBA_X(1) = -RNUM*((DOT_PROD(1,2)/DOT_PROD(1,1))*XR(1,1+1) - XR(1+1,1+2))/DEN

      FBA_Y(1) = -RNUM*((DOT_PROD(1,2)/DOT_PROD(1,1))*YR(1,1+1) - YR(1+1,1+2))/DEN

      FBA_Z(1) = -RNUM*((DOT_PROD(1,2)/DOT_PROD(1,1))*ZR(1,1+1) - ZR(1+1,1+2))/DEN

C }}}
C PARTICLE 2
C {{{
!     I = 2
      DEN = SINBOND(2)
      DEN1 = SINBOND(3)

      A1 = -RK_THETA*(BOND_ANGLE(2) - THETA_0)*( (DOT_PROD(2-1,2)/
     1  DOT_PROD(2,1))*XR(2,2+1) - (DOT_PROD(2-1,2)/DOT_PROD(2-1,1))
     1      *XR(2-1,2) + XR(2,2+1) - XR(2-1,2))/DEN

      A2 = -RK_THETA*(BOND_ANGLE(2+1) - THETA_0)*((DOT_PROD(2,2)/DOT_PROD(2,1))*XR(2,2+1) - XR(2+1,2+2))/DEN1

      FBA_X(2) = A1 + A2 

      A1 = -RK_THETA*(BOND_ANGLE(2) - THETA_0)*( (DOT_PROD(2-1,2)/
     1  DOT_PROD(2,1))*YR(2,2+1) - (DOT_PROD(2-1,2)/DOT_PROD(2-1,1))
     1      *YR(2-1,2) + YR(2,2+1) - YR(2-1,2))/DEN

      A2 = -RK_THETA*(BOND_ANGLE(2+1) - THETA_0)*((DOT_PROD(2,2)/DOT_PROD(2,1))*YR(2,2+1) - YR(2+1,2+2))/DEN1

      FBA_Y(2) = A1 + A2 

      A1 = -RK_THETA*(BOND_ANGLE(2) - THETA_0)*( (DOT_PROD(2-1,2)/
     1  DOT_PROD(2,1))*ZR(2,2+1) - (DOT_PROD(2-1,2)/DOT_PROD(2-1,1))
     1      *ZR(2-1,2) + ZR(2,2+1) - ZR(2-1,2))/DEN

      A2 = -RK_THETA*(BOND_ANGLE(2+1) - THETA_0)*((DOT_PROD(2,2)/DOT_PROD(2,1))*ZR(2,2+1) - ZR(2+1,2+2))/DEN1

      FBA_Z(2) = A1 + A2 

C }}}
C PARTICLES 3 THRU N-2 
C {{{
      DO I = 3, N-2

         DEN = SINBOND(I)
         DEN1 = SINBOND(I+1)
         DEN2 = SINBOND(I-1)

         A1 = -RK_THETA*(BOND_ANGLE(I) - THETA_0)*( (DOT_PROD(I-1,2)/
     1  DOT_PROD(I,1))*XR(I,I+1) - (DOT_PROD(I-1,2)/DOT_PROD(I-1,1))
     1      *XR(I-1,I) + XR(I,I+1) - XR(I-1,I))/DEN

         A2 = -RK_THETA*(BOND_ANGLE(I+1) - THETA_0)*((DOT_PROD(I,2)/DOT_PROD(I,1))*XR(I,I+1) - XR(I+1,I+2))/DEN1

         A3 = RK_THETA*(BOND_ANGLE(I-1) - THETA_0)*((DOT_PROD(I-2,2)/DOT_PROD(I-1,1))*XR(I-1,I) - XR(I-2,I-1))/DEN2

         FBA_X(I) = A1 + A2 + A3 

         A1 = -RK_THETA*(BOND_ANGLE(I) - THETA_0)*( (DOT_PROD(I-1,2)/
     1   DOT_PROD(I,1))*YR(I,I+1) - (DOT_PROD(I-1,2)/DOT_PROD(I-1,1))
     1      *YR(I-1,I) + YR(I,I+1) - YR(I-1,I))/DEN

         A2 = -RK_THETA*(BOND_ANGLE(I+1) - THETA_0)*((DOT_PROD(I,2)/
     1      DOT_PROD(I,1))*YR(I,I+1) - YR(I+1,I+2))/DEN1

         A3 = RK_THETA*(BOND_ANGLE(I-1) - THETA_0)*((DOT_PROD(I-2,2)/
     1      DOT_PROD(I-1,1))*YR(I-1,I) - YR(I-2,I-1))/DEN2

         FBA_Y(I) = A1 + A2 + A3 

         A1 = -RK_THETA*(BOND_ANGLE(I) - THETA_0)*( (DOT_PROD(I-1,2)/
     1   DOT_PROD(I,1))*ZR(I,I+1) - (DOT_PROD(I-1,2)/DOT_PROD(I-1,1))
     1      *ZR(I-1,I) + ZR(I,I+1) - ZR(I-1,I))/DEN

         A2 = -RK_THETA*(BOND_ANGLE(I+1) - THETA_0)*((DOT_PROD(I,2)/
     1      DOT_PROD(I,1))*ZR(I,I+1) - ZR(I+1,I+2))/DEN1

         A3 = RK_THETA*(BOND_ANGLE(I-1) - THETA_0)*((DOT_PROD(I-2,2)/
     1      DOT_PROD(I-1,1))*ZR(I-1,I) - ZR(I-2,I-1))/DEN2

         FBA_Z(I) = A1 + A2 + A3 

      ENDDO
C }}}
C PARTICLE N-1 
C {{{
!     I = N-1
      DEN = SINBOND(N-1)
      DEN1 = SINBOND(N-1-1)

      A1 = -RK_THETA*(BOND_ANGLE(N-1) - THETA_0)*( (DOT_PROD(N-1-1,2)/
     1   DOT_PROD(N-1,1))*XR(N-1,N-1+1) - (DOT_PROD(N-1-1,2)/DOT_PROD(N-1-1,1))
     1      *XR(N-1-1,N-1) + XR(N-1,N-1+1) - XR(N-1-1,N-1))/DEN

      A2 = RK_THETA*(BOND_ANGLE(N-1-1) - THETA_0)*((DOT_PROD(N-1-2,2)/
     1      DOT_PROD(N-1-1,1))*XR(N-1-1,N-1) - XR(N-1-2,N-1-1))/DEN1

      FBA_X(N-1) = A1 + A2

      A1 = -RK_THETA*(BOND_ANGLE(N-1) - THETA_0)*( (DOT_PROD(N-1-1,2)/
     1   DOT_PROD(N-1,1))*YR(N-1,N-1+1) - (DOT_PROD(N-1-1,2)/DOT_PROD(N-1-1,1))
     1      *YR(N-1-1,N-1) + YR(N-1,N-1+1) - YR(N-1-1,N-1))/DEN
   
      A2 = RK_THETA*(BOND_ANGLE(N-1-1) - THETA_0)*((DOT_PROD(N-1-2,2)/
     1      DOT_PROD(N-1-1,1))*YR(N-1-1,N-1) - YR(N-1-2,N-1-1))/DEN1

      FBA_Y(N-1) = A1 + A2

      A1 = -RK_THETA*(BOND_ANGLE(N-1) - THETA_0)*( (DOT_PROD(N-1-1,2)/
     1  DOT_PROD(N-1,1))*ZR(N-1,N-1+1) - (DOT_PROD(N-1-1,2)/DOT_PROD(N-1-1,1))
     1      *ZR(N-1-1,N-1) + ZR(N-1,N-1+1) - ZR(N-1-1,N-1))/DEN

      A2 = RK_THETA*(BOND_ANGLE(N-1-1) - THETA_0)*((DOT_PROD(N-1-2,2)/
     1      DOT_PROD(N-1-1,1))*ZR(N-1-1,N-1) - ZR(N-1-2,N-1-1))/DEN1

      FBA_Z(N-1) = A1 + A2
C }}}
C PARTICLE N
C {{{
!     I = N
      DEN = SINBOND(N-1)

      FBA_X(N) = RK_THETA*(BOND_ANGLE(N-1) - THETA_0)*
     1      ((DOT_PROD(N-2,2)/DOT_PROD(N-1,1))*XR(N-1,N) 
     1      - XR(N-2,N-1))/DEN

      FBA_Y(N) = RK_THETA*(BOND_ANGLE(N-1) - THETA_0)*
     1      ((DOT_PROD(N-2,2)/DOT_PROD(N-1,1))*YR(N-1,N) 
     1      - YR(N-2,N-1))/DEN

      FBA_Z(N) = RK_THETA*(BOND_ANGLE(N-1) - THETA_0)*
     1      ((DOT_PROD(N-2,2)/DOT_PROD(N-1,1))*ZR(N-1,N) 
     1      - ZR(N-2,N-1))/DEN

C }}}
C TORSIONAL ANGLE FORCES
C PARTICLES 1, 2, 3, N-2, N-1, AND N ARE DONE OUTSIDE OF THE LOOP
C PARTICLE 1
C {{{
!     I = 1
      COEF =DFAC(1+1)

      FTA_X(1) = -COEF*(-DOT_PROD(1+1,2)*XR(1+1,1+2) +
     1       DOT_PROD(1+1,1)*XR(1+2,1+3) -
     1      (1.0/X_PROD(1))*(DOT_PROD(1+1,2)*DOT_PROD(1,2) -
     1      DOT_PROD(1,3)*DOT_PROD(1+1,1))*(-DOT_PROD(1+1,1)*XR(1,1+1) +
     1      DOT_PROD(1,2)*XR(1+1,1+2))) 

      FTA_Y(1) = -COEF*(-DOT_PROD(1+1,2)*YR(1+1,1+2) +
     1      DOT_PROD(1+1,1)*YR(1+2,1+3) -
     1      (1.0/X_PROD(1))*(DOT_PROD(1+1,2)*DOT_PROD(1,2) -
     1      DOT_PROD(1,3)*DOT_PROD(1+1,1))*(-DOT_PROD(1+1,1)*YR(1,1+1) +
     1      DOT_PROD(1,2)*YR(1+1,1+2))) 

      FTA_Z(1) = -COEF*(-DOT_PROD(1+1,2)*ZR(1+1,1+2) +
     1      DOT_PROD(1+1,1)*ZR(1+2,1+3) -
     1      (1.0/X_PROD(1))*(DOT_PROD(1+1,2)*DOT_PROD(1,2) -
     1      DOT_PROD(1,3)*DOT_PROD(1+1,1))*(-DOT_PROD(1+1,1)*ZR(1,1+1) +
     1      DOT_PROD(1,2)*ZR(1+1,1+2))) 

C }}}
C PARTICLE 2
C {{{
!     I = 2
      COEF =DFAC(2+1)

      COEF1 = DFAC(2)

      A1 =  -COEF*(-DOT_PROD(2+1,2)*XR(2+1,2+2) +
     1      DOT_PROD(2+1,1)*XR(2+2,2+3) -
     1      (1.0/X_PROD(2))*(DOT_PROD(2+1,2)*DOT_PROD(2,2) -
     1      DOT_PROD(2,3)*DOT_PROD(2+1,1))*(-DOT_PROD(2+1,1)*XR(2,2+1) +
     1      DOT_PROD(2,2)*XR(2+1,2+2))) 

      A2 = -COEF1*(-DOT_PROD(2-1,2)*XR(2+1,2+2) +
     1      DOT_PROD(2,2)*XR(2,2+1) - DOT_PROD(2,2)*XR(2-1,2) -
     1      DOT_PROD(2,1)*XR(2+1,2+2) + 2.0*DOT_PROD(2-1,3)*XR(2,2+1) -
     1      (1.0/X_PROD(2-1))*(DOT_PROD(2,2)*DOT_PROD(2-1,2) -
     1      DOT_PROD(2-1,3)*DOT_PROD(2,1))*(DOT_PROD(2,1)*XR(2-1,2) -
     1      DOT_PROD(2-1,1)*XR(2,2+1) - DOT_PROD(2-1,2)*XR(2,2+1) +
     1      DOT_PROD(2-1,2)*XR(2-1,2)) -
     1      (1.0/X_PROD(2))*(DOT_PROD(2,2)*DOT_PROD(2-1,2) -
     1      DOT_PROD(2-1,3)*DOT_PROD(2,1))*(-DOT_PROD(2+1,1)*XR(2,2+1) +
     1      DOT_PROD(2,2)*XR(2+1,2+2))) 

      FTA_X(2) = A1 + A2 

      A1 = -DOT_PROD(3,2)*YR(3,4) + DOT_PROD(3,1)*YR(4,5) 
      A1=A1-(DOT_PROD(3,2)*DOT_PROD(2,2) - DOT_PROD(2,3)*DOT_PROD(3,1))*
     &     (-DOT_PROD(3,1)*YR(2,3) + DOT_PROD(2,2)*YR(3,4))/X_PROD(2)
      A1=-COEF*A1

      A2 = -COEF1*(-DOT_PROD(2-1,2)*YR(2+1,2+2) +
     1      DOT_PROD(2,2)*YR(2,2+1) - DOT_PROD(2,2)*YR(2-1,2) -
     1      DOT_PROD(2,1)*YR(2+1,2+2) + 2.0*DOT_PROD(2-1,3)*YR(2,2+1) -
     1      (1.0/X_PROD(2-1))*(DOT_PROD(2,2)*DOT_PROD(2-1,2) -
     1      DOT_PROD(2-1,3)*DOT_PROD(2,1))*(DOT_PROD(2,1)*YR(2-1,2) -
     1      DOT_PROD(2-1,1)*YR(2,2+1) - DOT_PROD(2-1,2)*YR(2,2+1) +
     1      DOT_PROD(2-1,2)*YR(2-1,2)) -
     1      (1.0/X_PROD(2))*(DOT_PROD(2,2)*DOT_PROD(2-1,2) -
     1      DOT_PROD(2-1,3)*DOT_PROD(2,1))*(-DOT_PROD(2+1,1)*YR(2,2+1) +
     1      DOT_PROD(2,2)*YR(2+1,2+2)))

      FTA_Y(2) = A1 + A2 
!     WRITE(MYUNIT,'(A,I6,5G15.5)') 'Y I,A1,A2,COEF,COEF1,FTA_Y=',2,A1,A2,COEF,COEF1,FTA_Y(2)
      
      A1 = -COEF*(-DOT_PROD(3,2)*ZR(3,4) +
     1      DOT_PROD(3,1)*ZR(4,5) -
     1      (1.0/X_PROD(2))*(DOT_PROD(3,2)*DOT_PROD(2,2) -
     1      DOT_PROD(2,3)*DOT_PROD(3,1))*(-DOT_PROD(3,1)*ZR(2,3) +
     1      DOT_PROD(2,2)*ZR(3,4))) 

      A2 = -COEF1*(-DOT_PROD(2-1,2)*ZR(2+1,2+2) +
     1      DOT_PROD(2,2)*ZR(2,2+1) - DOT_PROD(2,2)*ZR(2-1,2) -
     1      DOT_PROD(2,1)*ZR(2+1,2+2) + 2.0*DOT_PROD(2-1,3)*ZR(2,2+1) -
     1      (1.0/X_PROD(2-1))*(DOT_PROD(2,2)*DOT_PROD(2-1,2) -
     1      DOT_PROD(2-1,3)*DOT_PROD(2,1))*(DOT_PROD(2,1)*ZR(2-1,2) -
     1      DOT_PROD(2-1,1)*ZR(2,2+1) - DOT_PROD(2-1,2)*ZR(2,2+1) +
     1      DOT_PROD(2-1,2)*ZR(2-1,2)) -
     1      (1.0/X_PROD(2))*(DOT_PROD(2,2)*DOT_PROD(2-1,2) -
     1      DOT_PROD(2-1,3)*DOT_PROD(2,1))*(-DOT_PROD(2+1,1)*ZR(2,2+1) +
     1      DOT_PROD(2,2)*ZR(2+1,2+2))) 

      FTA_Z(2) = A1 + A2 
!     WRITE(MYUNIT,'(A,I6,4G15.5)') 'Z I,A1,A2,COEF,COEF1=',2,A1,A2,COEF,COEF1
C }}}
C PARTICLE 3
C {{{
!     I = 3
      COEF=DFAC(3+1)

      COEF1=DFAC(3)

      COEF2=DFAC(3-1)

      A1 = -COEF*(-DOT_PROD(3+1,2)*XR(3+1,3+2) +
     1      DOT_PROD(3+1,1)*XR(3+2,3+3) -
     1      (1.0/X_PROD(3))*(DOT_PROD(3+1,2)*DOT_PROD(3,2) -
     1      DOT_PROD(3,3)*DOT_PROD(3+1,1))*(-DOT_PROD(3+1,1)*XR(3,3+1) +
     1      DOT_PROD(3,2)*XR(3+1,3+2))) 

      A2 = -COEF1*(-DOT_PROD(3-1,2)*XR(3+1,3+2) +
     1      DOT_PROD(3,2)*XR(3,3+1) - DOT_PROD(3,2)*XR(3-1,3) -
     1      DOT_PROD(3,1)*XR(3+1,3+2) + 2.0*DOT_PROD(3-1,3)*XR(3,3+1) -
     1      (1.0/X_PROD(3-1))*(DOT_PROD(3,2)*DOT_PROD(3-1,2) -
     1      DOT_PROD(3-1,3)*DOT_PROD(3,1))*(DOT_PROD(3,1)*XR(3-1,3) -
     1      DOT_PROD(3-1,1)*XR(3,3+1) - DOT_PROD(3-1,2)*XR(3,3+1) +
     1      DOT_PROD(3-1,2)*XR(3-1,3)) -
     1      (1.0/X_PROD(3))*(DOT_PROD(3,2)*DOT_PROD(3-1,2) -
     1      DOT_PROD(3-1,3)*DOT_PROD(3,1))*(-DOT_PROD(3+1,1)*XR(3,3+1) +
     1      DOT_PROD(3,2)*XR(3+1,3+2))) 

      A3 = -COEF2*(DOT_PROD(3-2,2)*XR(3,3+1) -
     1      DOT_PROD(3-2,2)*XR(3-1,3) + DOT_PROD(3-1,2)*XR(3-2,3-1) +
     1      DOT_PROD(3-1,1)*XR(3-2,3-1) - 2.0*DOT_PROD(3-2,3)*XR(3-1,3) -
     1      (1.0/X_PROD(3-2))*(DOT_PROD(3-1,2)*DOT_PROD(3-2,2) -
     1      DOT_PROD(3-2,3)*DOT_PROD(3-1,1))*(DOT_PROD(3-2,1)*XR(3-1,3) -
     1      DOT_PROD(3-2,2)*XR(3-2,3-1)) -
     1      (1.0/X_PROD(3-1))*(DOT_PROD(3-1,2)*DOT_PROD(3-2,2) -
     1      DOT_PROD(3-2,3)*DOT_PROD(3-1,1))*(DOT_PROD(3,1)*XR(3-1,3) -
     1      DOT_PROD(3-1,1)*XR(3,3+1) - DOT_PROD(3-1,2)*XR(3,3+1) +
     1      DOT_PROD(3-1,2)*XR(3-1,3))) 

      FTA_X(3) = A1 + A2 + A3 
 
      A1 = -COEF*(-DOT_PROD(3+1,2)*YR(3+1,3+2) +
     1      DOT_PROD(3+1,1)*YR(3+2,3+3) -
     1      (1.0/X_PROD(3))*(DOT_PROD(3+1,2)*DOT_PROD(3,2) -
     1      DOT_PROD(3,3)*DOT_PROD(3+1,1))*(-DOT_PROD(3+1,1)*YR(3,3+1) +
     1      DOT_PROD(3,2)*YR(3+1,3+2))) 
      
      A2 = -COEF1*(-DOT_PROD(3-1,2)*YR(3+1,3+2) +
     1      DOT_PROD(3,2)*YR(3,3+1) - DOT_PROD(3,2)*YR(3-1,3) -
     1      DOT_PROD(3,1)*YR(3+1,3+2) + 2.0*DOT_PROD(3-1,3)*YR(3,3+1) -
     1      (1.0/X_PROD(3-1))*(DOT_PROD(3,2)*DOT_PROD(3-1,2) -
     1      DOT_PROD(3-1,3)*DOT_PROD(3,1))*(DOT_PROD(3,1)*YR(3-1,3) -
     1      DOT_PROD(3-1,1)*YR(3,3+1) - DOT_PROD(3-1,2)*YR(3,3+1) +
     1      DOT_PROD(3-1,2)*YR(3-1,3)) -
     1      (1.0/X_PROD(3))*(DOT_PROD(3,2)*DOT_PROD(3-1,2) -
     1      DOT_PROD(3-1,3)*DOT_PROD(3,1))*(-DOT_PROD(3+1,1)*YR(3,3+1) +
     1      DOT_PROD(3,2)*YR(3+1,3+2))) 

      A3 = -COEF2*(DOT_PROD(3-2,2)*YR(3,3+1) -
     1      DOT_PROD(3-2,2)*YR(3-1,3) + DOT_PROD(3-1,2)*YR(3-2,3-1) +
     1      DOT_PROD(3-1,1)*YR(3-2,3-1) - 2.0*DOT_PROD(3-2,3)*YR(3-1,3) -
     1      (1.0/X_PROD(3-2))*(DOT_PROD(3-1,2)*DOT_PROD(3-2,2) -
     1      DOT_PROD(3-2,3)*DOT_PROD(3-1,1))*(DOT_PROD(3-2,1)*YR(3-1,3) -
     1      DOT_PROD(3-2,2)*YR(3-2,3-1)) -
     1      (1.0/X_PROD(3-1))*(DOT_PROD(3-1,2)*DOT_PROD(3-2,2) -
     1      DOT_PROD(3-2,3)*DOT_PROD(3-1,1))*(DOT_PROD(3,1)*YR(3-1,3) -
     1      DOT_PROD(3-1,1)*YR(3,3+1) - DOT_PROD(3-1,2)*YR(3,3+1) +
     1      DOT_PROD(3-1,2)*YR(3-1,3)))

      FTA_Y(3) = A1 + A2 + A3 
 
      A1 = -COEF*(-DOT_PROD(3+1,2)*ZR(3+1,3+2) +
     1      DOT_PROD(3+1,1)*ZR(3+2,3+3) -
     1      (1.0/X_PROD(3))*(DOT_PROD(3+1,2)*DOT_PROD(3,2) -
     1      DOT_PROD(3,3)*DOT_PROD(3+1,1))*(-DOT_PROD(3+1,1)*ZR(3,3+1) +
     1      DOT_PROD(3,2)*ZR(3+1,3+2))) 

      A2 =  -COEF1*(-DOT_PROD(3-1,2)*ZR(3+1,3+2) +
     1      DOT_PROD(3,2)*ZR(3,3+1) - DOT_PROD(3,2)*ZR(3-1,3) -
     1      DOT_PROD(3,1)*ZR(3+1,3+2) + 2.0*DOT_PROD(3-1,3)*ZR(3,3+1) -
     1      (1.0/X_PROD(3-1))*(DOT_PROD(3,2)*DOT_PROD(3-1,2) -
     1      DOT_PROD(3-1,3)*DOT_PROD(3,1))*(DOT_PROD(3,1)*ZR(3-1,3) -
     1      DOT_PROD(3-1,1)*ZR(3,3+1) - DOT_PROD(3-1,2)*ZR(3,3+1) +
     1      DOT_PROD(3-1,2)*ZR(3-1,3)) -
     1      (1.0/X_PROD(3))*(DOT_PROD(3,2)*DOT_PROD(3-1,2) -
     1      DOT_PROD(3-1,3)*DOT_PROD(3,1))*(-DOT_PROD(3+1,1)*ZR(3,3+1) +
     1      DOT_PROD(3,2)*ZR(3+1,3+2))) 

      A3 = -COEF2*(DOT_PROD(3-2,2)*ZR(3,3+1) -
     1      DOT_PROD(3-2,2)*ZR(3-1,3) + DOT_PROD(3-1,2)*ZR(3-2,3-1) +
     1      DOT_PROD(3-1,1)*ZR(3-2,3-1) - 2.0*DOT_PROD(3-2,3)*ZR(3-1,3) -
     1      (1.0/X_PROD(3-2))*(DOT_PROD(3-1,2)*DOT_PROD(3-2,2) -
     1      DOT_PROD(3-2,3)*DOT_PROD(3-1,1))*(DOT_PROD(3-2,1)*ZR(3-1,3) -
     1      DOT_PROD(3-2,2)*ZR(3-2,3-1)) -
     1      (1.0/X_PROD(3-1))*(DOT_PROD(3-1,2)*DOT_PROD(3-2,2) -
     1      DOT_PROD(3-2,3)*DOT_PROD(3-1,1))*(DOT_PROD(3,1)*ZR(3-1,3) -
     1      DOT_PROD(3-1,1)*ZR(3,3+1) - DOT_PROD(3-1,2)*ZR(3,3+1) +
     1      DOT_PROD(3-1,2)*ZR(3-1,3))) 

      FTA_Z(3) = A1 + A2 + A3 
C }}}
C PARTICLES 4 TO N-3
C {{{
      DO I = 4, N-3

         COEF=DFAC(I+1)

         COEF1=DFAC(I)

         COEF2=DFAC(I-1)

         COEF3=DFAC(I-2)

         A1 = -COEF*(-DOT_PROD(I+1,2)*XR(I+1,I+2) +
     1      DOT_PROD(I+1,1)*XR(I+2,I+3) -
     1      (1.0/X_PROD(I))*(DOT_PROD(I+1,2)*DOT_PROD(I,2) -
     1      DOT_PROD(I,3)*DOT_PROD(I+1,1))*(-DOT_PROD(I+1,1)*XR(I,I+1) +
     1      DOT_PROD(I,2)*XR(I+1,I+2))) 

         A2 = -COEF1*(-DOT_PROD(I-1,2)*XR(I+1,I+2) +
     1      DOT_PROD(I,2)*XR(I,I+1) - DOT_PROD(I,2)*XR(I-1,I) -
     1      DOT_PROD(I,1)*XR(I+1,I+2) + 2.0*DOT_PROD(I-1,3)*XR(I,I+1) -
     1      (1.0/X_PROD(I-1))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1      DOT_PROD(I-1,3)*DOT_PROD(I,1))*(DOT_PROD(I,1)*XR(I-1,I) -
     1      DOT_PROD(I-1,1)*XR(I,I+1) - DOT_PROD(I-1,2)*XR(I,I+1) +
     1      DOT_PROD(I-1,2)*XR(I-1,I)) -
     1      (1.0/X_PROD(I))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1      DOT_PROD(I-1,3)*DOT_PROD(I,1))*(-DOT_PROD(I+1,1)*XR(I,I+1) +
     1      DOT_PROD(I,2)*XR(I+1,I+2))) 

         A3 = -COEF2*(DOT_PROD(I-2,2)*XR(I,I+1) -
     1      DOT_PROD(I-2,2)*XR(I-1,I) + DOT_PROD(I-1,2)*XR(I-2,I-1) +
     1      DOT_PROD(I-1,1)*XR(I-2,I-1) - 2.0*DOT_PROD(I-2,3)*XR(I-1,I) -
     1      (1.0/X_PROD(I-2))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1      DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I-2,1)*XR(I-1,I) -
     1      DOT_PROD(I-2,2)*XR(I-2,I-1)) -
     1      (1.0/X_PROD(I-1))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1  DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I,1)*XR(I-1,I) -
     1      DOT_PROD(I-1,1)*XR(I,I+1) - DOT_PROD(I-1,2)*XR(I,I+1) +
     1      DOT_PROD(I-1,2)*XR(I-1,I))) 

         A4 = -COEF3*(DOT_PROD(I-3,2)*XR(I-2,I-1) -
     1      DOT_PROD(I-2,1)*XR(I-3,I-2) -
     1      (1.0/X_PROD(I-2))*(DOT_PROD(I-2,2)*DOT_PROD(I-3,2) -
     1      DOT_PROD(I-3,3)*DOT_PROD(I-2,1))*(DOT_PROD(I-2,1)*XR(I-1,I) -
     1      DOT_PROD(I-2,2)*XR(I-2,I-1))) 

         FTA_X(I) = A1 + A2 + A3 + A4 

         A1 = -COEF*(-DOT_PROD(I+1,2)*YR(I+1,I+2) +
     1      DOT_PROD(I+1,1)*YR(I+2,I+3) -
     1      (1.0/X_PROD(I))*(DOT_PROD(I+1,2)*DOT_PROD(I,2) -
     1      DOT_PROD(I,3)*DOT_PROD(I+1,1))*(-DOT_PROD(I+1,1)*YR(I,I+1) +
     1      DOT_PROD(I,2)*YR(I+1,I+2))) 

         A2 = -COEF1*(-DOT_PROD(I-1,2)*YR(I+1,I+2) +
     1      DOT_PROD(I,2)*YR(I,I+1) - DOT_PROD(I,2)*YR(I-1,I) -
     1      DOT_PROD(I,1)*YR(I+1,I+2) + 2.0*DOT_PROD(I-1,3)*YR(I,I+1) -
     1      (1.0/X_PROD(I-1))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1      DOT_PROD(I-1,3)*DOT_PROD(I,1))*(DOT_PROD(I,1)*YR(I-1,I) -
     1      DOT_PROD(I-1,1)*YR(I,I+1) - DOT_PROD(I-1,2)*YR(I,I+1) +
     1      DOT_PROD(I-1,2)*YR(I-1,I)) -
     1      (1.0/X_PROD(I))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1      DOT_PROD(I-1,3)*DOT_PROD(I,1))*(-DOT_PROD(I+1,1)*YR(I,I+1) +
     1      DOT_PROD(I,2)*YR(I+1,I+2))) 

         A3 = -COEF2*(DOT_PROD(I-2,2)*YR(I,I+1) -
     1      DOT_PROD(I-2,2)*YR(I-1,I) + DOT_PROD(I-1,2)*YR(I-2,I-1) +
     1      DOT_PROD(I-1,1)*YR(I-2,I-1) - 2.0*DOT_PROD(I-2,3)*YR(I-1,I) -
     1      (1.0/X_PROD(I-2))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1      DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I-2,1)*YR(I-1,I) -
     1      DOT_PROD(I-2,2)*YR(I-2,I-1)) -
     1      (1.0/X_PROD(I-1))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1      DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I,1)*YR(I-1,I) -
     1      DOT_PROD(I-1,1)*YR(I,I+1) - DOT_PROD(I-1,2)*YR(I,I+1) +
     1      DOT_PROD(I-1,2)*YR(I-1,I))) 

         A4 = -COEF3*(DOT_PROD(I-3,2)*YR(I-2,I-1) -
     1      DOT_PROD(I-2,1)*YR(I-3,I-2) -
     1      (1.0/X_PROD(I-2))*(DOT_PROD(I-2,2)*DOT_PROD(I-3,2) -
     1      DOT_PROD(I-3,3)*DOT_PROD(I-2,1))*(DOT_PROD(I-2,1)*YR(I-1,I) -
     1      DOT_PROD(I-2,2)*YR(I-2,I-1))) 

         FTA_Y(I) = A1 + A2 + A3 + A4 

         A1 = -COEF*(-DOT_PROD(I+1,2)*ZR(I+1,I+2) +
     1      DOT_PROD(I+1,1)*ZR(I+2,I+3) -
     1      (1.0/X_PROD(I))*(DOT_PROD(I+1,2)*DOT_PROD(I,2) -
     1      DOT_PROD(I,3)*DOT_PROD(I+1,1))*(-DOT_PROD(I+1,1)*ZR(I,I+1) +
     1      DOT_PROD(I,2)*ZR(I+1,I+2))) 

         A2 = -COEF1*(-DOT_PROD(I-1,2)*ZR(I+1,I+2) +
     1      DOT_PROD(I,2)*ZR(I,I+1) - DOT_PROD(I,2)*ZR(I-1,I) -
     1      DOT_PROD(I,1)*ZR(I+1,I+2) + 2.0*DOT_PROD(I-1,3)*ZR(I,I+1) -
     1      (1.0/X_PROD(I-1))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1      DOT_PROD(I-1,3)*DOT_PROD(I,1))*(DOT_PROD(I,1)*ZR(I-1,I) -
     1      DOT_PROD(I-1,1)*ZR(I,I+1) - DOT_PROD(I-1,2)*ZR(I,I+1) +
     1      DOT_PROD(I-1,2)*ZR(I-1,I)) -
     1      (1.0/X_PROD(I))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1      DOT_PROD(I-1,3)*DOT_PROD(I,1))*(-DOT_PROD(I+1,1)*ZR(I,I+1) +
     1      DOT_PROD(I,2)*ZR(I+1,I+2))) 

         A3 = -COEF2*(DOT_PROD(I-2,2)*ZR(I,I+1) -
     1      DOT_PROD(I-2,2)*ZR(I-1,I) + DOT_PROD(I-1,2)*ZR(I-2,I-1) +
     1      DOT_PROD(I-1,1)*ZR(I-2,I-1) - 2.0*DOT_PROD(I-2,3)*ZR(I-1,I) -
     1      (1.0/X_PROD(I-2))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1      DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I-2,1)*ZR(I-1,I) -
     1      DOT_PROD(I-2,2)*ZR(I-2,I-1)) -
     1      (1.0/X_PROD(I-1))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1      DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I,1)*ZR(I-1,I) -
     1      DOT_PROD(I-1,1)*ZR(I,I+1) - DOT_PROD(I-1,2)*ZR(I,I+1) +
     1      DOT_PROD(I-1,2)*ZR(I-1,I))) 
      
         A4 = -COEF3*(DOT_PROD(I-3,2)*ZR(I-2,I-1) -
     1      DOT_PROD(I-2,1)*ZR(I-3,I-2) -
     1      (1.0/X_PROD(I-2))*(DOT_PROD(I-2,2)*DOT_PROD(I-3,2) -
     1      DOT_PROD(I-3,3)*DOT_PROD(I-2,1))*(DOT_PROD(I-2,1)*ZR(I-1,I) -
     1      DOT_PROD(I-2,2)*ZR(I-2,I-1))) 

         FTA_Z(I) = A1 + A2 + A3 + A4 

      ENDDO
C }}}
C PARTICLE N-2
C {{{
!     I = N-2
      COEF1=DFAC(N-2)

      COEF2=DFAC(N-2-1)

      COEF3=DFAC(N-2-2)

      A1 = -COEF1*(-DOT_PROD(N-2-1,2)*XR(N-2+1,N-2+2) + 
     1      DOT_PROD(N-2,2)*XR(N-2,N-2+1) - DOT_PROD(N-2,2)*XR(N-2-1,N-2) -
     1      DOT_PROD(N-2,1)*XR(N-2+1,N-2+2) + 2.0*DOT_PROD(N-2-1,3)*XR(N-2,N-2+1) -
     1      (1.0/X_PROD(N-2-1))*(DOT_PROD(N-2,2)*DOT_PROD(N-2-1,2) -
     1      DOT_PROD(N-2-1,3)*DOT_PROD(N-2,1))*(DOT_PROD(N-2,1)*XR(N-2-1,N-2) -
     1      DOT_PROD(N-2-1,1)*XR(N-2,N-2+1) - DOT_PROD(N-2-1,2)*XR(N-2,N-2+1) +
     1      DOT_PROD(N-2-1,2)*XR(N-2-1,N-2)) -
     1      (1.0/X_PROD(N-2))*(DOT_PROD(N-2,2)*DOT_PROD(N-2-1,2) -
     1      DOT_PROD(N-2-1,3)*DOT_PROD(N-2,1))*(-DOT_PROD(N-2+1,1)*XR(N-2,N-2+1) +
     1      DOT_PROD(N-2,2)*XR(N-2+1,N-2+2)))


      A2 = -COEF2*(DOT_PROD(N-2-2,2)*XR(N-2,N-2+1) -
     1      DOT_PROD(N-2-2,2)*XR(N-2-1,N-2) + DOT_PROD(N-2-1,2)*XR(N-2-2,N-2-1) +
     1      DOT_PROD(N-2-1,1)*XR(N-2-2,N-2-1) - 2.0*DOT_PROD(N-2-2,3)*XR(N-2-1,N-2) -
     1      (1.0/X_PROD(N-2-2))*(DOT_PROD(N-2-1,2)*DOT_PROD(N-2-2,2) -
     1      DOT_PROD(N-2-2,3)*DOT_PROD(N-2-1,1))*(DOT_PROD(N-2-2,1)*XR(N-2-1,N-2) -
     1      DOT_PROD(N-2-2,2)*XR(N-2-2,N-2-1)) -
     1      (1.0/X_PROD(N-2-1))*(DOT_PROD(N-2-1,2)*DOT_PROD(N-2-2,2) -
     1      DOT_PROD(N-2-2,3)*DOT_PROD(N-2-1,1))*(DOT_PROD(N-2,1)*XR(N-2-1,N-2) -
     1      DOT_PROD(N-2-1,1)*XR(N-2,N-2+1) - DOT_PROD(N-2-1,2)*XR(N-2,N-2+1) +
     1      DOT_PROD(N-2-1,2)*XR(N-2-1,N-2))) 
      
      A3 = -COEF3*(DOT_PROD(N-2-3,2)*XR(N-2-2,N-2-1) -
     1      DOT_PROD(N-2-2,1)*XR(N-2-3,N-2-2) -
     1      (1.0/X_PROD(N-2-2))*(DOT_PROD(N-2-2,2)*DOT_PROD(N-2-3,2) -
     1      DOT_PROD(N-2-3,3)*DOT_PROD(N-2-2,1))*(DOT_PROD(N-2-2,1)*XR(N-2-1,N-2) -
     1      DOT_PROD(N-2-2,2)*XR(N-2-2,N-2-1))) 

      FTA_X(N-2) = A1 + A2 + A3 

      A1 =  -COEF1*(-DOT_PROD(N-2-1,2)*YR(N-2+1,N-2+2) +  
     1      DOT_PROD(N-2,2)*YR(N-2,N-2+1) - DOT_PROD(N-2,2)*YR(N-2-1,N-2) -
     1      DOT_PROD(N-2,1)*YR(N-2+1,N-2+2) + 2.0*DOT_PROD(N-2-1,3)*YR(N-2,N-2+1) -
     1      (1.0/X_PROD(N-2-1))*(DOT_PROD(N-2,2)*DOT_PROD(N-2-1,2) -
     1      DOT_PROD(N-2-1,3)*DOT_PROD(N-2,1))*(DOT_PROD(N-2,1)*YR(N-2-1,N-2) -
     1      DOT_PROD(N-2-1,1)*YR(N-2,N-2+1) - DOT_PROD(N-2-1,2)*YR(N-2,N-2+1) +
     1      DOT_PROD(N-2-1,2)*YR(N-2-1,N-2)) -
     1      (1.0/X_PROD(N-2))*(DOT_PROD(N-2,2)*DOT_PROD(N-2-1,2) -
     1      DOT_PROD(N-2-1,3)*DOT_PROD(N-2,1))*(-DOT_PROD(N-2+1,1)*YR(N-2,N-2+1) +
     1      DOT_PROD(N-2,2)*YR(N-2+1,N-2+2))) 

      A2 =  -COEF2*(DOT_PROD(N-2-2,2)*YR(N-2,N-2+1) -
     1      DOT_PROD(N-2-2,2)*YR(N-2-1,N-2) + DOT_PROD(N-2-1,2)*YR(N-2-2,N-2-1) +
     1      DOT_PROD(N-2-1,1)*YR(N-2-2,N-2-1) - 2.0*DOT_PROD(N-2-2,3)*YR(N-2-1,N-2) -
     1      (1.0/X_PROD(N-2-2))*(DOT_PROD(N-2-1,2)*DOT_PROD(N-2-2,2) -
     1      DOT_PROD(N-2-2,3)*DOT_PROD(N-2-1,1))*(DOT_PROD(N-2-2,1)*YR(N-2-1,N-2) -
     1      DOT_PROD(N-2-2,2)*YR(N-2-2,N-2-1)) -
     1      (1.0/X_PROD(N-2-1))*(DOT_PROD(N-2-1,2)*DOT_PROD(N-2-2,2) -
     1      DOT_PROD(N-2-2,3)*DOT_PROD(N-2-1,1))*(DOT_PROD(N-2,1)*YR(N-2-1,N-2) -
     1      DOT_PROD(N-2-1,1)*YR(N-2,N-2+1) - DOT_PROD(N-2-1,2)*YR(N-2,N-2+1) +
     1      DOT_PROD(N-2-1,2)*YR(N-2-1,N-2)))

      A3 = -COEF3*(DOT_PROD(N-2-3,2)*YR(N-2-2,N-2-1) -
     1      DOT_PROD(N-2-2,1)*YR(N-2-3,N-2-2) -
     1      (1.0/X_PROD(N-2-2))*(DOT_PROD(N-2-2,2)*DOT_PROD(N-2-3,2) -
     1      DOT_PROD(N-2-3,3)*DOT_PROD(N-2-2,1))*(DOT_PROD(N-2-2,1)*YR(N-2-1,N-2) -
     1      DOT_PROD(N-2-2,2)*YR(N-2-2,N-2-1))) 

      FTA_Y(N-2) = A1 + A2 + A3 
 
      A1 = -COEF1*(-DOT_PROD(N-2-1,2)*ZR(N-2+1,N-2+2) +  
     1      DOT_PROD(N-2,2)*ZR(N-2,N-2+1) - DOT_PROD(N-2,2)*ZR(N-2-1,N-2) -
     1      DOT_PROD(N-2,1)*ZR(N-2+1,N-2+2) + 2.0*DOT_PROD(N-2-1,3)*ZR(N-2,N-2+1) -
     1      (1.0/X_PROD(N-2-1))*(DOT_PROD(N-2,2)*DOT_PROD(N-2-1,2) -
     1      DOT_PROD(N-2-1,3)*DOT_PROD(N-2,1))*(DOT_PROD(N-2,1)*ZR(N-2-1,N-2) -
     1      DOT_PROD(N-2-1,1)*ZR(N-2,N-2+1) - DOT_PROD(N-2-1,2)*ZR(N-2,N-2+1) +
     1      DOT_PROD(N-2-1,2)*ZR(N-2-1,N-2)) -
     1      (1.0/X_PROD(N-2))*(DOT_PROD(N-2,2)*DOT_PROD(N-2-1,2) -
     1      DOT_PROD(N-2-1,3)*DOT_PROD(N-2,1))*(-DOT_PROD(N-2+1,1)*ZR(N-2,N-2+1) +
     1      DOT_PROD(N-2,2)*ZR(N-2+1,N-2+2))) 

      A2 = -COEF2*(DOT_PROD(N-2-2,2)*ZR(N-2,N-2+1) -
     1      DOT_PROD(N-2-2,2)*ZR(N-2-1,N-2) + DOT_PROD(N-2-1,2)*ZR(N-2-2,N-2-1) +
     1      DOT_PROD(N-2-1,1)*ZR(N-2-2,N-2-1) - 2.0*DOT_PROD(N-2-2,3)*ZR(N-2-1,N-2) -
     1      (1.0/X_PROD(N-2-2))*(DOT_PROD(N-2-1,2)*DOT_PROD(N-2-2,2) -
     1      DOT_PROD(N-2-2,3)*DOT_PROD(N-2-1,1))*(DOT_PROD(N-2-2,1)*ZR(N-2-1,N-2) -
     1      DOT_PROD(N-2-2,2)*ZR(N-2-2,N-2-1)) -
     1      (1.0/X_PROD(N-2-1))*(DOT_PROD(N-2-1,2)*DOT_PROD(N-2-2,2) -
     1      DOT_PROD(N-2-2,3)*DOT_PROD(N-2-1,1))*(DOT_PROD(N-2,1)*ZR(N-2-1,N-2) -
     1      DOT_PROD(N-2-1,1)*ZR(N-2,N-2+1) - DOT_PROD(N-2-1,2)*ZR(N-2,N-2+1) +
     1      DOT_PROD(N-2-1,2)*ZR(N-2-1,N-2))) 

      A3 = -COEF3*(DOT_PROD(N-2-3,2)*ZR(N-2-2,N-2-1) -
     1      DOT_PROD(N-2-2,1)*ZR(N-2-3,N-2-2) -
     1      (1.0/X_PROD(N-2-2))*(DOT_PROD(N-2-2,2)*DOT_PROD(N-2-3,2) -
     1      DOT_PROD(N-2-3,3)*DOT_PROD(N-2-2,1))*(DOT_PROD(N-2-2,1)*ZR(N-2-1,N-2) -
     1      DOT_PROD(N-2-2,2)*ZR(N-2-2,N-2-1))) 

      FTA_Z(N-2) = A1 + A2 + A3 
C }}}
C PARTICLE N-1
C {{{
!     I = N-1
      COEF2=DFAC(N-1-1)

      COEF3=DFAC(N-1-2)

      A1 = -COEF2*(DOT_PROD(N-1-2,2)*XR(N-1,N-1+1) - 
     1      DOT_PROD(N-1-2,2)*XR(N-1-1,N-1) +
     1      DOT_PROD(N-1-1,2)*XR(N-1-2,N-1-1) +  DOT_PROD(N-1-1,1)*XR(N-1-2,N-1-1) -
     1      2.0*DOT_PROD(N-1-2,3)*XR(N-1-1,N-1) -
     1      (1.0/X_PROD(N-1-2))*(DOT_PROD(N-1-1,2)*DOT_PROD(N-1-2,2) -
     1      DOT_PROD(N-1-2,3)*DOT_PROD(N-1-1,1))*(DOT_PROD(N-1-2,1)*XR(N-1-1,N-1) -
     1      DOT_PROD(N-1-2,2)*XR(N-1-2,N-1-1)) -
     1      (1.0/X_PROD(N-1-1))*(DOT_PROD(N-1-1,2)*DOT_PROD(N-1-2,2) -
     1      DOT_PROD(N-1-2,3)*DOT_PROD(N-1-1,1))*(DOT_PROD(N-1,1)*XR(N-1-1,N-1) -
     1      DOT_PROD(N-1-1,1)*XR(N-1,N-1+1) - DOT_PROD(N-1-1,2)*XR(N-1,N-1+1) +
     1      DOT_PROD(N-1-1,2)*XR(N-1-1,N-1))) 

      A2 = -COEF3*(DOT_PROD(N-1-3,2)*XR(N-1-2,N-1-1) - 
     1      DOT_PROD(N-1-2,1)*XR(N-1-3,N-1-2) -
     1      (1.0/X_PROD(N-1-2))*(DOT_PROD(N-1-2,2)*DOT_PROD(N-1-3,2) -
     1      DOT_PROD(N-1-3,3)*DOT_PROD(N-1-2,1))*(DOT_PROD(N-1-2,1)*XR(N-1-1,N-1) -
     1      DOT_PROD(N-1-2,2)*XR(N-1-2,N-1-1))) 

      FTA_X(N-1) = A1 + A2  

      A1 = -COEF2*(DOT_PROD(N-1-2,2)*YR(N-1,N-1+1) - 
     1      DOT_PROD(N-1-2,2)*YR(N-1-1,N-1) +
     1      DOT_PROD(N-1-1,2)*YR(N-1-2,N-1-1) +  DOT_PROD(N-1-1,1)*YR(N-1-2,N-1-1) -
     1      2.0*DOT_PROD(N-1-2,3)*YR(N-1-1,N-1) -
     1      (1.0/X_PROD(N-1-2))*(DOT_PROD(N-1-1,2)*DOT_PROD(N-1-2,2) -
     1      DOT_PROD(N-1-2,3)*DOT_PROD(N-1-1,1))*(DOT_PROD(N-1-2,1)*YR(N-1-1,N-1) -
     1      DOT_PROD(N-1-2,2)*YR(N-1-2,N-1-1)) -
     1      (1.0/X_PROD(N-1-1))*(DOT_PROD(N-1-1,2)*DOT_PROD(N-1-2,2) -
     1      DOT_PROD(N-1-2,3)*DOT_PROD(N-1-1,1))*(DOT_PROD(N-1,1)*YR(N-1-1,N-1) -
     1  DOT_PROD(N-1-1,1)*YR(N-1,N-1+1) - DOT_PROD(N-1-1,2)*YR(N-1,N-1+1) +
     1      DOT_PROD(N-1-1,2)*YR(N-1-1,N-1))) 

      A2 = -COEF3*(DOT_PROD(N-1-3,2)*YR(N-1-2,N-1-1) - 
     1      DOT_PROD(N-1-2,1)*YR(N-1-3,N-1-2) -
     1      (1.0/X_PROD(N-1-2))*(DOT_PROD(N-1-2,2)*DOT_PROD(N-1-3,2) -
     1      DOT_PROD(N-1-3,3)*DOT_PROD(N-1-2,1))*(DOT_PROD(N-1-2,1)*YR(N-1-1,N-1) -
     1      DOT_PROD(N-1-2,2)*YR(N-1-2,N-1-1))) 

      FTA_Y(N-1) = A1 + A2  

      A1 = -COEF2*(DOT_PROD(N-1-2,2)*ZR(N-1,N-1+1) - 
     1      DOT_PROD(N-1-2,2)*ZR(N-1-1,N-1) +
     1      DOT_PROD(N-1-1,2)*ZR(N-1-2,N-1-1) +  DOT_PROD(N-1-1,1)*ZR(N-1-2,N-1-1) -
     1      2.0*DOT_PROD(N-1-2,3)*ZR(N-1-1,N-1) -
     1      (1.0/X_PROD(N-1-2))*(DOT_PROD(N-1-1,2)*DOT_PROD(N-1-2,2) -
     1      DOT_PROD(N-1-2,3)*DOT_PROD(N-1-1,1))*(DOT_PROD(N-1-2,1)*ZR(N-1-1,N-1) -
     1      DOT_PROD(N-1-2,2)*ZR(N-1-2,N-1-1)) -
     1      (1.0/X_PROD(N-1-1))*(DOT_PROD(N-1-1,2)*DOT_PROD(N-1-2,2) -
     1      DOT_PROD(N-1-2,3)*DOT_PROD(N-1-1,1))*(DOT_PROD(N-1,1)*ZR(N-1-1,N-1) -
     1      DOT_PROD(N-1-1,1)*ZR(N-1,N-1+1) - DOT_PROD(N-1-1,2)*ZR(N-1,N-1+1) +
     1      DOT_PROD(N-1-1,2)*ZR(N-1-1,N-1))) 

      A2 = -COEF3*(DOT_PROD(N-1-3,2)*ZR(N-1-2,N-1-1) - 
     1      DOT_PROD(N-1-2,1)*ZR(N-1-3,N-1-2) -
     1      (1.0/X_PROD(N-1-2))*(DOT_PROD(N-1-2,2)*DOT_PROD(N-1-3,2) -
     1      DOT_PROD(N-1-3,3)*DOT_PROD(N-1-2,1))*(DOT_PROD(N-1-2,1)*ZR(N-1-1,N-1) -
     1      DOT_PROD(N-1-2,2)*ZR(N-1-2,N-1-1))) 

      FTA_Z(N-1) = A1 + A2 
C }}} 
C PARTICLE N
C {{{
!     I = N
      COEF3=DFAC(N-2)

      FTA_X(N) = -COEF3*(DOT_PROD(N-3,2)*XR(N-2,N-1) 
     1      - DOT_PROD(N-2,1)*XR(N-3,N-2) -
     1      (1.0/X_PROD(N-2))*(DOT_PROD(N-2,2)*DOT_PROD(N-3,2) -
     1      DOT_PROD(N-3,3)*DOT_PROD(N-2,1))*(DOT_PROD(N-2,1)*XR(N-1,N) -
     1      DOT_PROD(N-2,2)*XR(N-2,N-1))) 

      FTA_Y(N) = -COEF3*(DOT_PROD(N-3,2)*YR(N-2,N-1) - 
     1      DOT_PROD(N-2,1)*YR(N-3,N-2) 
     1      - (1.0/X_PROD(N-2))*(DOT_PROD(N-2,2)*DOT_PROD(N-3,2) 
     1      - DOT_PROD(N-3,3)*DOT_PROD(N-2,1))*(DOT_PROD(N-2,1)*YR(N-1,N) -
     1      DOT_PROD(N-2,2)*YR(N-2,N-1))) 

      FTA_Z(N) = -COEF3*(DOT_PROD(N-3,2)*ZR(N-2,N-1) - 
     1      DOT_PROD(N-2,1)*ZR(N-3,N-2) -
     1      (1.0/X_PROD(N-2))*(DOT_PROD(N-2,2)*DOT_PROD(N-3,2) -
     1      DOT_PROD(N-3,3)*DOT_PROD(N-2,1))*(DOT_PROD(N-2,1)*ZR(N-1,N) -
     1      DOT_PROD(N-2,2)*ZR(N-2,N-1))) 
C }}}
C TOTAL UP THE GRADIENTS
C {{{
      DO I = 1, N
!        IF (I.EQ.2) THEN
!           WRITE(MYUNIT,'(A,I6,4F15.5)') 'I,FNBX,FBX,FBAX,FTAX=',I,FNB_X(I),FB_X(I),FBA_X(I),FTA_X(I)
!           WRITE(MYUNIT,'(A,I6,4F15.5)') 'I,FNBY,FBY,FBAY,FTAY=',I,FNB_Y(I),FB_Y(I),FBA_Y(I),FTA_Y(I)
!           WRITE(MYUNIT,'(A,I6,4F15.5)') 'I,FNBZ,FBZ,FBAZ,FTAZ=',I,FNB_Z(I),FB_Z(I),FBA_Z(I),FTA_Z(I)
!        ENDIF
         FX(I) = FNB_X(I) + FB_X(I) + FBA_X(I) + FTA_X(I) 
         FY(I) = FNB_Y(I) + FB_Y(I) + FBA_Y(I) + FTA_Y(I) 
         FZ(I) = FNB_Z(I) + FB_Z(I) + FBA_Z(I) + FTA_Z(I) 
      ENDDO

      DO I = 1, N
         J = (I-1)*3
         FQ(J+1) = -FX(I)
         FQ(J+2) = -FY(I)
         FQ(J+3) = -FZ(I)
      ENDDO
C }}}
      RETURN
      END
C }}}
C
C>  FILL THE PARAMETER ARRAYS
C
      SUBROUTINE PARAM_ARRAYBLN(LJREP_BLN,LJATT_BLN,A_BLN,B_BLN,C_BLN,D_BLN,BEADLETTER,BLNSSTRUCT,
     &   LJREPBB, LJATTBB, LJREPLL, LJATTLL, LJREPNN, LJATTNN, 
     &   HABLN, HBBLN, HCBLN, HDBLN, EABLN, EBBLN, ECBLN, EDBLN, TABLN, TBBLN, TCBLN, TDBLN, N)
C {{{
C DECLARATIONS {{{
      IMPLICIT NONE
      INTEGER N
      INTEGER NTYPE(N), J1, I, J, ICOUNT
      DOUBLE PRECISION LJREP_BLN(N,N), LJATT_BLN(N,N), A_BLN(N), B_BLN(N), C_BLN(N), D_BLN(N),
     &                 LJREPBB, LJATTBB, LJREPLL, LJATTLL, LJREPNN, LJATTNN, 
     &                 HABLN, HBBLN, HCBLN, HDBLN, EABLN, EBBLN, ECBLN, EDBLN, TABLN, TBBLN, TCBLN, TDBLN
      CHARACTER(LEN=1) BEADLETTER(N), BLNSSTRUCT(N)
C }}}
C AMINO ACID TYPES {{{
C B=1 L=2 N=3
C
      DO J1=1,N
         IF (BEADLETTER(J1).EQ.'B') THEN
            NTYPE(J1)=1
         ELSEIF (BEADLETTER(J1).EQ.'L') THEN
            NTYPE(J1)=2
         ELSEIF (BEADLETTER(J1).EQ.'N') THEN
            NTYPE(J1)=3
         ELSE
            PRINT '(A,A1)','ERROR IN PARAM_ARRAYBLN, UNRECOGNISED BEAD TYPE: ',BEADLETTER(J1)
            STOP
         ENDIF
      ENDDO

C }}}

C PARAMETERS FOR THE DIHEDRAL ANGLE POTENTIAL 
C THE END BONDS HAVE NO DIHEDRAL TERM, SO THE TOTAL NUMBER OF TERMS
C IS N-3 (E.G. 4 ATOMS, 1 DIHEDRAL). NON-ZERO TERMS FOR
C 2 TO 3, 3 TO 4, 4 TO 5, ... , N-3 TO N-2, N-2 TO N-1. THE
C H, E AND T PARAMETERS ARE DEFINED FOR THE FIRST BEAD OF EACH EDGE,
C I.E. FOR 2, 3, 4, ..., N-2.
C {{{

      A_BLN(1:N)=0.0D0
      B_BLN(1:N)=0.0D0
      C_BLN(1:N)=0.0D0
      D_BLN(1:N)=0.0D0
      DO I=1,N-3
         IF (BLNSSTRUCT(I).EQ.'H') THEN
            A_BLN(I+1)=HABLN
            B_BLN(I+1)=HBBLN
            C_BLN(I+1)=HCBLN
            D_BLN(I+1)=HDBLN
         ELSE IF (BLNSSTRUCT(I).EQ.'E') THEN
            A_BLN(I+1)=EABLN
            B_BLN(I+1)=EBBLN
            C_BLN(I+1)=ECBLN
            D_BLN(I+1)=EDBLN
         ELSE IF (BLNSSTRUCT(I).EQ.'T') THEN
            A_BLN(I+1)=TABLN
            B_BLN(I+1)=TBBLN
            C_BLN(I+1)=TCBLN
            D_BLN(I+1)=TDBLN
         ELSE
            PRINT '(A,A1)','ERROR IN PARAM_ARRAYBLN, UNRECOGNISED SS TYPE: ',BLNSSTRUCT(J1)
            STOP
         ENDIF
C        PRINT '(A,I6,A,A1,A,4F12.4)','I+1=',I+1,' SYMBOL=',BLNSSTRUCT(I),' A,B,C,D=',A_BLN(I+1),B_BLN(I+1),C_BLN(I+1),D_BLN(I+1)
      ENDDO
C }}}
C  PARAMETERS FOR THE L-J INTERACTION BETWEEN NON-BONDED PARTICLES {{{

      DO I = 1, N-1
         DO J = I+1, N
            IF (NTYPE(I) .EQ. 3 .OR. NTYPE(J) .EQ. 3) THEN
               LJREP_BLN(I,J) = LJREPNN
               LJATT_BLN(I,J) = LJATTNN
               LJREP_BLN(J,I) = LJREPNN
               LJATT_BLN(J,I) = LJATTNN
            ELSE IF (NTYPE(I) .EQ. 1 .AND. NTYPE(J) .EQ. 1) THEN
               LJREP_BLN(I,J) = LJREPBB
               LJATT_BLN(I,J) = LJATTBB
               LJREP_BLN(J,I) = LJREPBB
               LJATT_BLN(J,I) = LJATTBB
            ELSE
               LJREP_BLN(I,J) = LJREPLL
               LJATT_BLN(I,J) = LJATTLL
               LJREP_BLN(J,I) = LJREPLL
               LJATT_BLN(J,I) = LJATTLL
            ENDIF
         ENDDO
      ENDDO
C }}}

      RETURN
      END
C }}}
