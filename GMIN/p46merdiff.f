C GPL LICENSE INFO{{{
C   OPTIM: A PROGRAM FOR OPTIMIZING GEOMETRIES AND CALCULATING REACTION PATHWAYS
C   COPYRIGHT (C) 1999-2006 DAVID J. WALES
C   THIS FILE IS PART OF OPTIM.
C
C   OPTIM IS FREE SOFTWARE; YOU CAN REDISTRIBUTE IT AND/OR MODIFY
C   IT UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE AS PUBLISHED BY
C   THE FREE SOFTWARE FOUNDATION; EITHER VERSION 2 OF THE LICENSE, OR
C   (AT YOUR OPTION) ANY LATER VERSION.
C
C   OPTIM IS DISTRIBUTED IN THE HOPE THAT IT WILL BE USEFUL,
C   BUT WITHOUT ANY WARRANTY; WITHOUT EVEN THE IMPLIED WARRANTY OF
C   MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  SEE THE
C   GNU GENERAL PUBLIC LICENSE FOR MORE DETAILS.
C
C   YOU SHOULD HAVE RECEIVED A COPY OF THE GNU GENERAL PUBLIC LICENSE
C   ALONG WITH THIS PROGRAM; IF NOT, WRITE TO THE FREE SOFTWARE
C   FOUNDATION, INC., 59 TEMPLE PLACE, SUITE 330, BOSTON, MA  02111-1307  USA
C
C}}}

C \MAINPAGE PROGRAM: P46MERDIFF.F
C>
C> CALCULATE THE ENERGY, GRADIENT, AND SECOND
C> DERIVATIVES FOR A GIVEN CONFIGURATION OF THE 46 PARTICLE POLYMER CHAIN.
C> A CONFIGURATION AND NUMBER OF PARTICLES IS PASSED TO THE SUBROUTINE AND
C> THE ENERGY (ENERGY), GRADIENT (GRAD), AND 
C> THE MATRIX OF SECOND DERIVATIVES ARE RETURNED.
C> AUTHOR: JOHN ROSE

        SUBROUTINE P46MERDIFF(QO, N, GRAD, ENERGY, GTEST)
C {{{
        IMPLICIT NONE
        INTEGER NTYPE(46),N
        DOUBLE PRECISION RMASS, EPSILON,SIGMA,DELTA,THETA_0,RK_R,RK_THETA,QO(3*N),GRAD(3*N),ENERGY
        PARAMETER (RMASS = 40.0, EPSILON = 0.0100570)
        PARAMETER (SIGMA=3.4, DELTA=1.0D-6, THETA_0 = 1.8326)
        PARAMETER (RK_R = 20.0*0.0100570, RK_THETA = 20.0*0.0100570)
        LOGICAL GTEST, STEST
        DOUBLE PRECISION A_PARAM(N,N), B_PARAM(N,N), D_PARAM(N),C_PARAM(N),
     1                  X(N), Y(N), Z(N), XR(N,N), YR(N,N), ZR(N,N),
     2                  DOT_PROD(N,3), X_PROD(N), BOND_ANGLE(N), TOR_ANGLE(N), RADII(N,N)

C       COMMON/WORK/A_PARAM(N,N),
C    1  B_PARAM(N,N),NTYPE(46),
C    2  D_PARAM(N),C_PARAM(N),
C    3  X(N), Y(N), Z(N), 
C    4  XR(N,N), YR(N,N), ZR(N,N), 
C    5  DOT_PROD(N,3), X_PROD(N), 
C    6  BOND_ANGLE(N), STEST, TOR_ANGLE(N), RADII(N,N)

        STEST=.FALSE.

        CALL PARAM_ARRAY(A_PARAM,B_PARAM,C_PARAM,D_PARAM,N)
        CALL CALC_INT_COORDS(QO,N,A_PARAM,B_PARAM,C_PARAM,D_PARAM,X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD,
     1                       BOND_ANGLE,TOR_ANGLE,RADII,NTYPE)
        CALL CALC_ENERGY(QO,ENERGY,N,A_PARAM,B_PARAM,C_PARAM,D_PARAM,X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD,
     1                   BOND_ANGLE,TOR_ANGLE,RADII,NTYPE)
        IF ((.NOT.GTEST).AND.(.NOT.STEST)) RETURN
        CALL CALC_GRADIENT(QO,GRAD,N,A_PARAM,B_PARAM,C_PARAM,D_PARAM,X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD,
     1                     BOND_ANGLE,TOR_ANGLE,RADII,NTYPE)

        IF (.NOT.STEST) RETURN
        CALL CALC_DYN(QO,N,A_PARAM,B_PARAM,C_PARAM,D_PARAM,X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD,
     1                BOND_ANGLE,TOR_ANGLE,RADII,NTYPE)

        RETURN
        END
C }}}
C> CALCULATE THE INTERNAL COORDINATES

        SUBROUTINE CALC_INT_COORDS(QO,N,A_PARAM,B_PARAM,C_PARAM,D_PARAM,X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD,
     1                             BOND_ANGLE,TOR_ANGLE,RADII,NTYPE)
C {{{
        IMPLICIT NONE
        INTEGER NTYPE(46),I,J,N
        DOUBLE PRECISION QO(3*N)
        DOUBLE PRECISION A_PARAM(N,N), B_PARAM(N,N), D_PARAM(N),C_PARAM(N),
     1       X(N), Y(N), Z(N), XR(N,N), YR(N,N), ZR(N,N),
     2       DOT_PROD(N,3), X_PROD(N), BOND_ANGLE(N), TOR_ANGLE(N), RADII(N,N),
     3       COS_THETA, COS_PHI

C       COMMON/WORK/A_PARAM(N,N),
C    1  B_PARAM(N,N),NTYPE(46),
C    2  D_PARAM(N),C_PARAM(N),
C    3  X(N), Y(N), Z(N),
C    4  XR(N,N), YR(N,N), ZR(N,N),
C    5  DOT_PROD(N,3), X_PROD(N),
C    6  BOND_ANGLE(N), TOR_ANGLE(N), RADII(N,N)

        DO I = 1, N
        J = (I-1)*3
        X(I) = QO(J+1)
        Y(I) = QO(J+2)
        Z(I) = QO(J+3)
        ENDDO

C INTER-PARTICLE DISTANCES

        DO I = 1, N-1
        DO J = I+1, N
        XR(I,J) = X(J) - X(I)
        YR(I,J) = Y(J) - Y(I)
        ZR(I,J) = Z(J) - Z(I)
        RADII(I,J) = DSQRT(XR(I,J)*XR(I,J) + YR(I,J)*YR(I,J) 
     1  + ZR(I,J)*ZR(I,J))
        RADII(J,I) = RADII(I,J)
        ENDDO
        ENDDO

C DOT PRODUCTS BETWEEN BOND VECTORS

        DO I = 1, N-3
        DOT_PROD(I,1) = XR(I,I+1)*XR(I,I+1) + YR(I,I+1)*YR(I,I+1) +
     1  ZR(I,I+1)*ZR(I,I+1)

        DOT_PROD(I,2) = XR(I,I+1)*XR(I+1,I+2)+YR(I,I+1)*YR(I+1,I+2)+ 
     1  ZR(I,I+1)*ZR(I+1,I+2)

        DOT_PROD(I,3) = XR(I,I+1)*XR(I+2,I+3)+YR(I,I+1)*YR(I+2,I+3)+
     1  ZR(I,I+1)*ZR(I+2,I+3)
        ENDDO

        I = N-2
        DOT_PROD(I,1) = XR(I,I+1)*XR(I,I+1) + YR(I,I+1)*YR(I,I+1) +
     1  ZR(I,I+1)*ZR(I,I+1)

        DOT_PROD(I,2) = XR(I,I+1)*XR(I+1,I+2)+YR(I,I+1)*YR(I+1,I+2)+
     1  ZR(I,I+1)*ZR(I+1,I+2)

        I = N-1
        DOT_PROD(I,1) = XR(I,I+1)*XR(I,I+1) + YR(I,I+1)*YR(I,I+1) +
     1  ZR(I,I+1)*ZR(I,I+1)

C CROSS-PRODUCTS BETWEEN ADJACENT BOND VECTORS

        DO I = 1, N-2
        X_PROD(I) = DOT_PROD(I,1)*DOT_PROD(I+1,1) -
     1               DOT_PROD(I,2)*DOT_PROD(I,2)   
        ENDDO

C BOND ANGLES

        DO I = 1, N-2
        COS_THETA=-DOT_PROD(I,2)/(DSQRT(DOT_PROD(I,1)
     1  *DOT_PROD(I+1,1)))
        BOND_ANGLE(I+1) = DACOS(COS_THETA)
        ENDDO

C TORSIONAL ANGLES

        DO I = 1, N-3
        COS_PHI = (DOT_PROD(I,2)*DOT_PROD(I+1,2) -
     1  DOT_PROD(I,3)*DOT_PROD(I+1,1))/DSQRT(X_PROD(I)*X_PROD(I+1))
        IF (ABS(COS_PHI).GT.1.0D0) COS_PHI=COS_PHI/ABS(COS_PHI)
        TOR_ANGLE(I+1) = DACOS(COS_PHI)
C       WRITE(*,'(A,I4,4F20.10)') 'I,TOR_ANGLE,COS_PHI,DACOS=',I,TOR_ANGLE(I+1),COS_PHI,DACOS(COS_PHI)
        ENDDO

        RETURN
        END

C }}} 
C> CALCULATE THE ENERGY
        SUBROUTINE CALC_ENERGY(QO,ENERGY,N,A_PARAM,B_PARAM,C_PARAM,D_PARAM,X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD,
     1                         BOND_ANGLE,TOR_ANGLE,RADII,NTYPE)
C {{{
        IMPLICIT NONE
        INTEGER I,J,N
        DOUBLE PRECISION RMASS, EPSILON, SIGMA, DELTA, THETA_0, RK_R, RK_THETA, RAD6, E_TANGLE,
     1                   QO(*), ENERGY, S6, E_NBOND, E_BOND, E_BANGLE
        PARAMETER (RMASS = 40.0, EPSILON = 0.0100570)
        PARAMETER (SIGMA=3.4, DELTA=1.0D-6, THETA_0 = 1.8326)
        PARAMETER (RK_R = 20.0*0.0100570, RK_THETA = 20.0*0.0100570)
        INTEGER NTYPE(46)
        DOUBLE PRECISION A_PARAM(N,N), B_PARAM(N,N), D_PARAM(N),C_PARAM(N),
     1                  X(N), Y(N), Z(N), XR(N,N), YR(N,N), ZR(N,N),
     2                  DOT_PROD(N,3), X_PROD(N), BOND_ANGLE(N), TOR_ANGLE(N), RADII(N,N)

C       COMMON/WORK/A_PARAM(N,N),
C    1  B_PARAM(N,N),NTYPE(46),
C    2  D_PARAM(N),C_PARAM(N),
C    3  X(N), Y(N), Z(N),
C    4  XR(N,N), YR(N,N), ZR(N,N),
C    5  DOT_PROD(N,3), X_PROD(N),
C    6  BOND_ANGLE(N), TOR_ANGLE(N), RADII(N,N)

        S6 = SIGMA*SIGMA*SIGMA*SIGMA*SIGMA*SIGMA
        E_NBOND=0.0D0
        E_BOND=0.0D0
        E_BANGLE=0.0D0
        E_TANGLE=0.0D0

        DO I = 1, N-2
        DO J = I+2, N

        RAD6 = RADII(I,J)*RADII(I,J)*RADII(I,J)*RADII(I,J)*RADII(I,J)*
     1  RADII(I,J)

        E_NBOND = E_NBOND + 4.0*((A_PARAM(I,J)*S6*S6/(RAD6*RAD6)) + 
     1  (B_PARAM(I,J)*S6/RAD6))

        ENDDO
        ENDDO

        DO I = 1, N-1

        E_BOND = E_BOND + 0.5*RK_R*(RADII(I,I+1)-SIGMA)*
     1  (RADII(I,I+1)-SIGMA)

        ENDDO

        
        DO I = 2, N-1

        E_BANGLE = E_BANGLE + 0.5*RK_THETA*(BOND_ANGLE(I)-THETA_0)
     1  *(BOND_ANGLE(I)-THETA_0)

        ENDDO

        DO I = 2, N-2

        E_TANGLE = E_TANGLE + C_PARAM(I)*(1.0 + COS(TOR_ANGLE(I))) 
     1  + D_PARAM(I)*(1.0 + COS(3.0*TOR_ANGLE(I)))

        ENDDO

        ENERGY = E_NBOND + E_BOND + E_BANGLE + E_TANGLE
C       WRITE(*,'(A,4F20.10)') 'NBOND,BOND,BANGLE,TANGLE=',E_NBOND,E_BOND,E_BANGLE,E_TANGLE

        RETURN
        END
C }}}
C CALCULATE THE GRADIENTS

        SUBROUTINE CALC_GRADIENT(QO,FQ,N,A_PARAM,B_PARAM,C_PARAM,D_PARAM,X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD, 
     1                           BOND_ANGLE,TOR_ANGLE,RADII,NTYPE)
C {{{
C DECLARATIONS {{{
        IMPLICIT NONE
        INTEGER I, J, N
        DOUBLE PRECISION RMASS, EPSILON, SIGMA, DELTA, RK_R, RK_THETA,
     1                   A4, COEF, COEF1, COEF2, COEF3, A3, DEN2, A2, A1, DEN1, RNUM, DEN,
     2                   RVAR, FZZ, FYY, FXX, DF, RAD14, RAD7, S6, QO(3*N), THETA_0
        PARAMETER (RMASS = 40.0, EPSILON = 0.0100570)
        PARAMETER (SIGMA=3.4, DELTA=1.0D-6,THETA_0 = 1.8326)
        PARAMETER (RK_R = 20.0*0.0100570, RK_THETA = 20.0*0.0100570)
        INTEGER NTYPE(46)
        DOUBLE PRECISION A_PARAM(N,N), B_PARAM(N,N), D_PARAM(N),C_PARAM(N),
     1                  X(N),Y(N),Z(N),XR(N,N), YR(N,N), ZR(N,N),
     2     DOT_PROD(N,3), X_PROD(N),BOND_ANGLE(N),TOR_ANGLE(N), RADII(N,N)
        DOUBLE PRECISION FNB_X(N),FNB_Y(N),FNB_Z(N),
     1  FB_X(N),FB_Y(N),FB_Z(N)
        DOUBLE PRECISION FBA_X(N),FBA_Y(N),FBA_Z(N), FX(N),FY(N),FZ(N)
        DOUBLE PRECISION FTA_X(N),FTA_Y(N),FTA_Z(N), FQ(3*N)
C }}}
C       COMMON/WORK/A_PARAM(N,N),
C    1  B_PARAM(N,N),NTYPE(46),
C    2  D_PARAM(N),C_PARAM(N),
C    3  X(N), Y(N), Z(N),
C    4  XR(N,N), YR(N,N), ZR(N,N),
C    5  DOT_PROD(N,3), X_PROD(N),
C    6  BOND_ANGLE(N), TOR_ANGLE(N), RADII(N,N)

        S6 = SIGMA*SIGMA*SIGMA*SIGMA*SIGMA*SIGMA

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
C ..... NON-BONDED INTERACTION FORCES ..... {{{

        DO I = 1, N-2
        DO J = I+2, N

        RAD7 = RADII(I,J)*RADII(I,J)*RADII(I,J)*RADII(I,J)*
     1        RADII(I,J)*RADII(I,J)*RADII(I,J)   
        RAD14 = RAD7*RAD7 

        DF = -24.0*((2.0*A_PARAM(I,J)*S6*S6/RAD14) + 
     1              (B_PARAM(I,J)*S6/(RAD7*RADII(I,J))))

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
C ... BOND INTERACTION FORCES ... {{{

        DO I = 1, N-1

        RVAR = SIGMA/RADII(I,I+1) 

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
        I = 1
        DEN = DSIN(BOND_ANGLE(I+1))
     1        *DSQRT(DOT_PROD(I+1,1)*DOT_PROD(I,1))
        RNUM = RK_THETA*(BOND_ANGLE(I+1) - THETA_0)

        FBA_X(I) = -RNUM*((DOT_PROD(I,2)/DOT_PROD(I,1))*XR(I,I+1) -
     1        XR(I+1,I+2))/DEN

        FBA_Y(I) = -RNUM*((DOT_PROD(I,2)/DOT_PROD(I,1))*YR(I,I+1) -
     1        YR(I+1,I+2))/DEN

        FBA_Z(I) = -RNUM*((DOT_PROD(I,2)/DOT_PROD(I,1))*ZR(I,I+1) -
     1        ZR(I+1,I+2))/DEN

C PARTICLE 2

        I = 2
        DEN = DSIN(BOND_ANGLE(I))
     1        *DSQRT(DOT_PROD(I,1)*DOT_PROD(I-1,1))
        DEN1 = DSIN(BOND_ANGLE(I+1))*DSQRT(DOT_PROD(I+1,1)
     1         *DOT_PROD(I,1))

        A1 = -RK_THETA*(BOND_ANGLE(I) - THETA_0)*( (DOT_PROD(I-1,2)/
     1  DOT_PROD(I,1))*XR(I,I+1) - (DOT_PROD(I-1,2)/DOT_PROD(I-1,1))
     1        *XR(I-1,I) + XR(I,I+1) - XR(I-1,I))/DEN

        A2 = -RK_THETA*(BOND_ANGLE(I+1) - THETA_0)*((DOT_PROD(I,2)/
     1        DOT_PROD(I,1))*XR(I,I+1) - XR(I+1,I+2))/DEN1

        FBA_X(I) = A1 + A2 

        A1 = -RK_THETA*(BOND_ANGLE(I) - THETA_0)*( (DOT_PROD(I-1,2)/
     1  DOT_PROD(I,1))*YR(I,I+1) - (DOT_PROD(I-1,2)/DOT_PROD(I-1,1))
     1        *YR(I-1,I) + YR(I,I+1) - YR(I-1,I))/DEN

        A2 = -RK_THETA*(BOND_ANGLE(I+1) - THETA_0)*((DOT_PROD(I,2)/
     1        DOT_PROD(I,1))*YR(I,I+1) - YR(I+1,I+2))/DEN1

        FBA_Y(I) = A1 + A2 

        A1 = -RK_THETA*(BOND_ANGLE(I) - THETA_0)*( (DOT_PROD(I-1,2)/
     1  DOT_PROD(I,1))*ZR(I,I+1) - (DOT_PROD(I-1,2)/DOT_PROD(I-1,1))
     1        *ZR(I-1,I) + ZR(I,I+1) - ZR(I-1,I))/DEN

        A2 = -RK_THETA*(BOND_ANGLE(I+1) - THETA_0)*((DOT_PROD(I,2)/
     1        DOT_PROD(I,1))*ZR(I,I+1) - ZR(I+1,I+2))/DEN1

        FBA_Z(I) = A1 + A2 

C PARTICLES 3 THRU N-2 

        DO I = 3, N-2

        DEN = DSIN(BOND_ANGLE(I))*
     1              DSQRT(DOT_PROD(I,1)*DOT_PROD(I-1,1))
        DEN1 = DSIN(BOND_ANGLE(I+1))*
     1               DSQRT(DOT_PROD(I+1,1)*DOT_PROD(I,1))
        DEN2 = DSIN(BOND_ANGLE(I-1))*DSQRT(DOT_PROD(I-2,1)
     1         *DOT_PROD(I-1,1))

        A1 = -RK_THETA*(BOND_ANGLE(I) - THETA_0)*( (DOT_PROD(I-1,2)/
     1  DOT_PROD(I,1))*XR(I,I+1) - (DOT_PROD(I-1,2)/DOT_PROD(I-1,1))
     1        *XR(I-1,I) + XR(I,I+1) - XR(I-1,I))/DEN

        A2 = -RK_THETA*(BOND_ANGLE(I+1) - THETA_0)*((DOT_PROD(I,2)/
     1        DOT_PROD(I,1))*XR(I,I+1) - XR(I+1,I+2))/DEN1

        A3 = RK_THETA*(BOND_ANGLE(I-1) - THETA_0)*((DOT_PROD(I-2,2)/
     1        DOT_PROD(I-1,1))*XR(I-1,I) - XR(I-2,I-1))/DEN2

        FBA_X(I) = A1 + A2 + A3 

        A1 = -RK_THETA*(BOND_ANGLE(I) - THETA_0)*( (DOT_PROD(I-1,2)/
     1  DOT_PROD(I,1))*YR(I,I+1) - (DOT_PROD(I-1,2)/DOT_PROD(I-1,1))
     1        *YR(I-1,I) + YR(I,I+1) - YR(I-1,I))/DEN

        A2 = -RK_THETA*(BOND_ANGLE(I+1) - THETA_0)*((DOT_PROD(I,2)/
     1        DOT_PROD(I,1))*YR(I,I+1) - YR(I+1,I+2))/DEN1

        A3 = RK_THETA*(BOND_ANGLE(I-1) - THETA_0)*((DOT_PROD(I-2,2)/
     1        DOT_PROD(I-1,1))*YR(I-1,I) - YR(I-2,I-1))/DEN2

        FBA_Y(I) = A1 + A2 + A3 

        A1 = -RK_THETA*(BOND_ANGLE(I) - THETA_0)*( (DOT_PROD(I-1,2)/
     1  DOT_PROD(I,1))*ZR(I,I+1) - (DOT_PROD(I-1,2)/DOT_PROD(I-1,1))
     1        *ZR(I-1,I) + ZR(I,I+1) - ZR(I-1,I))/DEN

        A2 = -RK_THETA*(BOND_ANGLE(I+1) - THETA_0)*((DOT_PROD(I,2)/
     1        DOT_PROD(I,1))*ZR(I,I+1) - ZR(I+1,I+2))/DEN1

        A3 = RK_THETA*(BOND_ANGLE(I-1) - THETA_0)*((DOT_PROD(I-2,2)/
     1        DOT_PROD(I-1,1))*ZR(I-1,I) - ZR(I-2,I-1))/DEN2

        FBA_Z(I) = A1 + A2 + A3 

        ENDDO

C PARTICLE N-1 

        I = N-1
        DEN = DSIN(BOND_ANGLE(I))*
     1              DSQRT(DOT_PROD(I,1)*DOT_PROD(I-1,1))
        DEN1 = DSIN(BOND_ANGLE(I-1))*DSQRT(DOT_PROD(I-2,1)
     1         *DOT_PROD(I-1,1))

        A1 = -RK_THETA*(BOND_ANGLE(I) - THETA_0)*( (DOT_PROD(I-1,2)/
     1  DOT_PROD(I,1))*XR(I,I+1) - (DOT_PROD(I-1,2)/DOT_PROD(I-1,1))
     1        *XR(I-1,I) + XR(I,I+1) - XR(I-1,I))/DEN

        A2 = RK_THETA*(BOND_ANGLE(I-1) - THETA_0)*((DOT_PROD(I-2,2)/
     1        DOT_PROD(I-1,1))*XR(I-1,I) - XR(I-2,I-1))/DEN1

        FBA_X(I) = A1 + A2

        A1 = -RK_THETA*(BOND_ANGLE(I) - THETA_0)*( (DOT_PROD(I-1,2)/
     1  DOT_PROD(I,1))*YR(I,I+1) - (DOT_PROD(I-1,2)/DOT_PROD(I-1,1))
     1        *YR(I-1,I) + YR(I,I+1) - YR(I-1,I))/DEN

        A2 = RK_THETA*(BOND_ANGLE(I-1) - THETA_0)*((DOT_PROD(I-2,2)/
     1        DOT_PROD(I-1,1))*YR(I-1,I) - YR(I-2,I-1))/DEN1

        FBA_Y(I) = A1 + A2

        A1 = -RK_THETA*(BOND_ANGLE(I) - THETA_0)*( (DOT_PROD(I-1,2)/
     1  DOT_PROD(I,1))*ZR(I,I+1) - (DOT_PROD(I-1,2)/DOT_PROD(I-1,1))
     1        *ZR(I-1,I) + ZR(I,I+1) - ZR(I-1,I))/DEN

        A2 = RK_THETA*(BOND_ANGLE(I-1) - THETA_0)*((DOT_PROD(I-2,2)/
     1        DOT_PROD(I-1,1))*ZR(I-1,I) - ZR(I-2,I-1))/DEN1

        FBA_Z(I) = A1 + A2

C PARTICLE N

        I = N
        DEN = DSIN(BOND_ANGLE(I-1))*DSQRT(DOT_PROD(I-2,1)
     1        *DOT_PROD(I-1,1))

        FBA_X(I) = RK_THETA*(BOND_ANGLE(I-1) - THETA_0)*
     1        ((DOT_PROD(I-2,2)/DOT_PROD(I-1,1))*XR(I-1,I) 
     1        - XR(I-2,I-1))/DEN

        FBA_Y(I) = RK_THETA*(BOND_ANGLE(I-1) - THETA_0)*
     1        ((DOT_PROD(I-2,2)/DOT_PROD(I-1,1))*YR(I-1,I) 
     1        - YR(I-2,I-1))/DEN

        FBA_Z(I) = RK_THETA*(BOND_ANGLE(I-1) - THETA_0)*
     1        ((DOT_PROD(I-2,2)/DOT_PROD(I-1,1))*ZR(I-1,I) 
     1        - ZR(I-2,I-1))/DEN
C }}}
C TORSIONAL ANGLE FORCES
C PARTICLES 1, 2, 3, N-2, N-1, AND N ARE DONE OUTSIDE OF THE LOOP
C PARTICLE 1

        I = 1
             COEF =(C_PARAM(I+1)+D_PARAM(I+1)*(12.0*DCOS(TOR_ANGLE(I+1))
     1         *DCOS(TOR_ANGLE(I+1))-3.0))
     1  *(1.0D0/DSQRT(X_PROD(I+1)*X_PROD(I)))  

        FTA_X(I) = -COEF*(-DOT_PROD(I+1,2)*XR(I+1,I+2) +
     1         DOT_PROD(I+1,1)*XR(I+2,I+3) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I+1,2)*DOT_PROD(I,2) -
     1        DOT_PROD(I,3)*DOT_PROD(I+1,1))*(-DOT_PROD(I+1,1)*XR(I,I+1) +
     1        DOT_PROD(I,2)*XR(I+1,I+2))) 

        FTA_Y(I) = -COEF*(-DOT_PROD(I+1,2)*YR(I+1,I+2) +
     1        DOT_PROD(I+1,1)*YR(I+2,I+3) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I+1,2)*DOT_PROD(I,2) -
     1        DOT_PROD(I,3)*DOT_PROD(I+1,1))*(-DOT_PROD(I+1,1)*YR(I,I+1) +
     1        DOT_PROD(I,2)*YR(I+1,I+2))) 

        FTA_Z(I) = -COEF*(-DOT_PROD(I+1,2)*ZR(I+1,I+2) +
     1        DOT_PROD(I+1,1)*ZR(I+2,I+3) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I+1,2)*DOT_PROD(I,2) -
     1        DOT_PROD(I,3)*DOT_PROD(I+1,1))*(-DOT_PROD(I+1,1)*ZR(I,I+1) +
     1        DOT_PROD(I,2)*ZR(I+1,I+2))) 


C PARTICLE 2

        I = 2
        COEF =(C_PARAM(I+1)+D_PARAM(I+1)*(12.0*DCOS(TOR_ANGLE(I+1))
     1  *DCOS(TOR_ANGLE(I+1)) - 3.0))
     1        *(1.0D0/DSQRT(X_PROD(I+1)*X_PROD(I)))  

             COEF1 = (C_PARAM(I) + D_PARAM(I)*(12.0*DCOS(TOR_ANGLE(I))
     1        *DCOS(TOR_ANGLE(I)) - 
     1        3.0))*(1.0D0/DSQRT(X_PROD(I)*X_PROD(I-1)))  

        A1 =  -COEF*(-DOT_PROD(I+1,2)*XR(I+1,I+2) +
     1        DOT_PROD(I+1,1)*XR(I+2,I+3) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I+1,2)*DOT_PROD(I,2) -
     1        DOT_PROD(I,3)*DOT_PROD(I+1,1))*(-DOT_PROD(I+1,1)*XR(I,I+1) +
     1        DOT_PROD(I,2)*XR(I+1,I+2))) 

        A2 = -COEF1*(-DOT_PROD(I-1,2)*XR(I+1,I+2) +
     1        DOT_PROD(I,2)*XR(I,I+1) - DOT_PROD(I,2)*XR(I-1,I) -
     1        DOT_PROD(I,1)*XR(I+1,I+2) + 2.0*DOT_PROD(I-1,3)*XR(I,I+1) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(DOT_PROD(I,1)*XR(I-1,I) -
     1        DOT_PROD(I-1,1)*XR(I,I+1) - DOT_PROD(I-1,2)*XR(I,I+1) +
     1        DOT_PROD(I-1,2)*XR(I-1,I)) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(-DOT_PROD(I+1,1)*XR(I,I+1) +
     1        DOT_PROD(I,2)*XR(I+1,I+2))) 

        FTA_X(I) = A1 + A2 

        A1 = -COEF*(-DOT_PROD(I+1,2)*YR(I+1,I+2) +
     1        DOT_PROD(I+1,1)*YR(I+2,I+3) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I+1,2)*DOT_PROD(I,2) -
     1        DOT_PROD(I,3)*DOT_PROD(I+1,1))*(-DOT_PROD(I+1,1)*YR(I,I+1) +
     1        DOT_PROD(I,2)*YR(I+1,I+2))) 

        A2 = -COEF1*(-DOT_PROD(I-1,2)*YR(I+1,I+2) +
     1        DOT_PROD(I,2)*YR(I,I+1) - DOT_PROD(I,2)*YR(I-1,I) -
     1        DOT_PROD(I,1)*YR(I+1,I+2) + 2.0*DOT_PROD(I-1,3)*YR(I,I+1) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(DOT_PROD(I,1)*YR(I-1,I) -
     1        DOT_PROD(I-1,1)*YR(I,I+1) - DOT_PROD(I-1,2)*YR(I,I+1) +
     1        DOT_PROD(I-1,2)*YR(I-1,I)) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(-DOT_PROD(I+1,1)*YR(I,I+1) +
     1        DOT_PROD(I,2)*YR(I+1,I+2)))

        FTA_Y(I) = A1 + A2 
        
        A1 = -COEF*(-DOT_PROD(I+1,2)*ZR(I+1,I+2) +
     1        DOT_PROD(I+1,1)*ZR(I+2,I+3) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I+1,2)*DOT_PROD(I,2) -
     1        DOT_PROD(I,3)*DOT_PROD(I+1,1))*(-DOT_PROD(I+1,1)*ZR(I,I+1) +
     1        DOT_PROD(I,2)*ZR(I+1,I+2))) 

        A2 = -COEF1*(-DOT_PROD(I-1,2)*ZR(I+1,I+2) +
     1        DOT_PROD(I,2)*ZR(I,I+1) - DOT_PROD(I,2)*ZR(I-1,I) -
     1        DOT_PROD(I,1)*ZR(I+1,I+2) + 2.0*DOT_PROD(I-1,3)*ZR(I,I+1) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(DOT_PROD(I,1)*ZR(I-1,I) -
     1        DOT_PROD(I-1,1)*ZR(I,I+1) - DOT_PROD(I-1,2)*ZR(I,I+1) +
     1        DOT_PROD(I-1,2)*ZR(I-1,I)) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(-DOT_PROD(I+1,1)*ZR(I,I+1) +
     1        DOT_PROD(I,2)*ZR(I+1,I+2))) 

        FTA_Z(I) = A1 + A2 

C PARTICLE 3

        I = 3
        COEF=(C_PARAM(I+1)+D_PARAM(I+1)*(12.0*DCOS(TOR_ANGLE(I+1))
     1  *DCOS(TOR_ANGLE(I+1)) - 
     1        3.0))*(1.0D0/DSQRT(X_PROD(I+1)*X_PROD(I)))  

        COEF1=(C_PARAM(I)+D_PARAM(I)*(12.0*DCOS(TOR_ANGLE(I))
     1        *DCOS(TOR_ANGLE(I)) - 
     1        3.0))*(1.0D0/DSQRT(X_PROD(I)*X_PROD(I-1)))  

        COEF2=(C_PARAM(I-1)+D_PARAM(I-1)*(12.0*DCOS(TOR_ANGLE(I-1))
     1        *DCOS(TOR_ANGLE(I-1)) - 
     1        3.0))*(1.0D0/DSQRT(X_PROD(I-1)*X_PROD(I-2)))  

        A1 = -COEF*(-DOT_PROD(I+1,2)*XR(I+1,I+2) +
     1        DOT_PROD(I+1,1)*XR(I+2,I+3) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I+1,2)*DOT_PROD(I,2) -
     1        DOT_PROD(I,3)*DOT_PROD(I+1,1))*(-DOT_PROD(I+1,1)*XR(I,I+1) +
     1        DOT_PROD(I,2)*XR(I+1,I+2))) 

        A2 = -COEF1*(-DOT_PROD(I-1,2)*XR(I+1,I+2) +
     1        DOT_PROD(I,2)*XR(I,I+1) - DOT_PROD(I,2)*XR(I-1,I) -
     1        DOT_PROD(I,1)*XR(I+1,I+2) + 2.0*DOT_PROD(I-1,3)*XR(I,I+1) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(DOT_PROD(I,1)*XR(I-1,I) -
     1        DOT_PROD(I-1,1)*XR(I,I+1) - DOT_PROD(I-1,2)*XR(I,I+1) +
     1        DOT_PROD(I-1,2)*XR(I-1,I)) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(-DOT_PROD(I+1,1)*XR(I,I+1) +
     1        DOT_PROD(I,2)*XR(I+1,I+2))) 

        A3 = -COEF2*(DOT_PROD(I-2,2)*XR(I,I+1) -
     1        DOT_PROD(I-2,2)*XR(I-1,I) + DOT_PROD(I-1,2)*XR(I-2,I-1) +
     1        DOT_PROD(I-1,1)*XR(I-2,I-1) - 2.0*DOT_PROD(I-2,3)*XR(I-1,I) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I-2,1)*XR(I-1,I) -
     1        DOT_PROD(I-2,2)*XR(I-2,I-1)) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I,1)*XR(I-1,I) -
     1        DOT_PROD(I-1,1)*XR(I,I+1) - DOT_PROD(I-1,2)*XR(I,I+1) +
     1        DOT_PROD(I-1,2)*XR(I-1,I))) 

        FTA_X(I) = A1 + A2 + A3 
 
        A1 = -COEF*(-DOT_PROD(I+1,2)*YR(I+1,I+2) +
     1        DOT_PROD(I+1,1)*YR(I+2,I+3) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I+1,2)*DOT_PROD(I,2) -
     1        DOT_PROD(I,3)*DOT_PROD(I+1,1))*(-DOT_PROD(I+1,1)*YR(I,I+1) +
     1        DOT_PROD(I,2)*YR(I+1,I+2))) 
        
        A2 = -COEF1*(-DOT_PROD(I-1,2)*YR(I+1,I+2) +
     1        DOT_PROD(I,2)*YR(I,I+1) - DOT_PROD(I,2)*YR(I-1,I) -
     1        DOT_PROD(I,1)*YR(I+1,I+2) + 2.0*DOT_PROD(I-1,3)*YR(I,I+1) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(DOT_PROD(I,1)*YR(I-1,I) -
     1        DOT_PROD(I-1,1)*YR(I,I+1) - DOT_PROD(I-1,2)*YR(I,I+1) +
     1        DOT_PROD(I-1,2)*YR(I-1,I)) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(-DOT_PROD(I+1,1)*YR(I,I+1) +
     1        DOT_PROD(I,2)*YR(I+1,I+2))) 

        A3 = -COEF2*(DOT_PROD(I-2,2)*YR(I,I+1) -
     1        DOT_PROD(I-2,2)*YR(I-1,I) + DOT_PROD(I-1,2)*YR(I-2,I-1) +
     1        DOT_PROD(I-1,1)*YR(I-2,I-1) - 2.0*DOT_PROD(I-2,3)*YR(I-1,I) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I-2,1)*YR(I-1,I) -
     1        DOT_PROD(I-2,2)*YR(I-2,I-1)) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I,1)*YR(I-1,I) -
     1        DOT_PROD(I-1,1)*YR(I,I+1) - DOT_PROD(I-1,2)*YR(I,I+1) +
     1        DOT_PROD(I-1,2)*YR(I-1,I)))

        FTA_Y(I) = A1 + A2 + A3 
 
        A1 = -COEF*(-DOT_PROD(I+1,2)*ZR(I+1,I+2) +
     1        DOT_PROD(I+1,1)*ZR(I+2,I+3) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I+1,2)*DOT_PROD(I,2) -
     1        DOT_PROD(I,3)*DOT_PROD(I+1,1))*(-DOT_PROD(I+1,1)*ZR(I,I+1) +
     1        DOT_PROD(I,2)*ZR(I+1,I+2))) 

        A2 =  -COEF1*(-DOT_PROD(I-1,2)*ZR(I+1,I+2) +
     1        DOT_PROD(I,2)*ZR(I,I+1) - DOT_PROD(I,2)*ZR(I-1,I) -
     1        DOT_PROD(I,1)*ZR(I+1,I+2) + 2.0*DOT_PROD(I-1,3)*ZR(I,I+1) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(DOT_PROD(I,1)*ZR(I-1,I) -
     1        DOT_PROD(I-1,1)*ZR(I,I+1) - DOT_PROD(I-1,2)*ZR(I,I+1) +
     1        DOT_PROD(I-1,2)*ZR(I-1,I)) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(-DOT_PROD(I+1,1)*ZR(I,I+1) +
     1        DOT_PROD(I,2)*ZR(I+1,I+2))) 

        A3 = -COEF2*(DOT_PROD(I-2,2)*ZR(I,I+1) -
     1        DOT_PROD(I-2,2)*ZR(I-1,I) + DOT_PROD(I-1,2)*ZR(I-2,I-1) +
     1        DOT_PROD(I-1,1)*ZR(I-2,I-1) - 2.0*DOT_PROD(I-2,3)*ZR(I-1,I) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I-2,1)*ZR(I-1,I) -
     1        DOT_PROD(I-2,2)*ZR(I-2,I-1)) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I,1)*ZR(I-1,I) -
     1        DOT_PROD(I-1,1)*ZR(I,I+1) - DOT_PROD(I-1,2)*ZR(I,I+1) +
     1        DOT_PROD(I-1,2)*ZR(I-1,I))) 

        FTA_Z(I) = A1 + A2 + A3 

C PARTICLES 4 TO N-3

        DO I = 4, N-3

        COEF = (C_PARAM(I+1) + D_PARAM(I+1)*(12.0*DCOS(TOR_ANGLE(I+1))
     1  *DCOS(TOR_ANGLE(I+1)) - 3.0))*(1.0D0/DSQRT(X_PROD(I+1)*X_PROD(I)))

        COEF1 = (C_PARAM(I) + D_PARAM(I)*(12.0*DCOS(TOR_ANGLE(I))
     1        *DCOS(TOR_ANGLE(I)) -3.0))*(1.0D0/DSQRT(X_PROD(I)*X_PROD(I-1)))  

        COEF2 = (C_PARAM(I-1) + D_PARAM(I-1)*(12.0*DCOS(TOR_ANGLE(I-1))
     1        *DCOS(TOR_ANGLE(I-1)) - 
     1  3.0))*(1.0D0/DSQRT(X_PROD(I-1)*X_PROD(I-2)))  

        COEF3 = (C_PARAM(I-2) + D_PARAM(I-2)*(12.0*DCOS(TOR_ANGLE(I-2))
     1        *DCOS(TOR_ANGLE(I-2)) - 
     1  3.0))*(1.0D0/DSQRT(X_PROD(I-2)*X_PROD(I-3)))  

        A1 = -COEF*(-DOT_PROD(I+1,2)*XR(I+1,I+2) +
     1        DOT_PROD(I+1,1)*XR(I+2,I+3) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I+1,2)*DOT_PROD(I,2) -
     1        DOT_PROD(I,3)*DOT_PROD(I+1,1))*(-DOT_PROD(I+1,1)*XR(I,I+1) +
     1        DOT_PROD(I,2)*XR(I+1,I+2))) 

        A2 = -COEF1*(-DOT_PROD(I-1,2)*XR(I+1,I+2) +
     1        DOT_PROD(I,2)*XR(I,I+1) - DOT_PROD(I,2)*XR(I-1,I) -
     1        DOT_PROD(I,1)*XR(I+1,I+2) + 2.0*DOT_PROD(I-1,3)*XR(I,I+1) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(DOT_PROD(I,1)*XR(I-1,I) -
     1        DOT_PROD(I-1,1)*XR(I,I+1) - DOT_PROD(I-1,2)*XR(I,I+1) +
     1        DOT_PROD(I-1,2)*XR(I-1,I)) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(-DOT_PROD(I+1,1)*XR(I,I+1) +
     1        DOT_PROD(I,2)*XR(I+1,I+2))) 

        A3 = -COEF2*(DOT_PROD(I-2,2)*XR(I,I+1) -
     1        DOT_PROD(I-2,2)*XR(I-1,I) + DOT_PROD(I-1,2)*XR(I-2,I-1) +
     1        DOT_PROD(I-1,1)*XR(I-2,I-1) - 2.0*DOT_PROD(I-2,3)*XR(I-1,I) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I-2,1)*XR(I-1,I) -
     1        DOT_PROD(I-2,2)*XR(I-2,I-1)) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1  DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I,1)*XR(I-1,I) -
     1        DOT_PROD(I-1,1)*XR(I,I+1) - DOT_PROD(I-1,2)*XR(I,I+1) +
     1        DOT_PROD(I-1,2)*XR(I-1,I))) 

        A4 = -COEF3*(DOT_PROD(I-3,2)*XR(I-2,I-1) -
     1        DOT_PROD(I-2,1)*XR(I-3,I-2) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-2,2)*DOT_PROD(I-3,2) -
     1        DOT_PROD(I-3,3)*DOT_PROD(I-2,1))*(DOT_PROD(I-2,1)*XR(I-1,I) -
     1        DOT_PROD(I-2,2)*XR(I-2,I-1))) 

        FTA_X(I) = A1 + A2 + A3 + A4 

        A1 = -COEF*(-DOT_PROD(I+1,2)*YR(I+1,I+2) +
     1        DOT_PROD(I+1,1)*YR(I+2,I+3) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I+1,2)*DOT_PROD(I,2) -
     1        DOT_PROD(I,3)*DOT_PROD(I+1,1))*(-DOT_PROD(I+1,1)*YR(I,I+1) +
     1        DOT_PROD(I,2)*YR(I+1,I+2))) 

        A2 = -COEF1*(-DOT_PROD(I-1,2)*YR(I+1,I+2) +
     1        DOT_PROD(I,2)*YR(I,I+1) - DOT_PROD(I,2)*YR(I-1,I) -
     1        DOT_PROD(I,1)*YR(I+1,I+2) + 2.0*DOT_PROD(I-1,3)*YR(I,I+1) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(DOT_PROD(I,1)*YR(I-1,I) -
     1        DOT_PROD(I-1,1)*YR(I,I+1) - DOT_PROD(I-1,2)*YR(I,I+1) +
     1        DOT_PROD(I-1,2)*YR(I-1,I)) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(-DOT_PROD(I+1,1)*YR(I,I+1) +
     1        DOT_PROD(I,2)*YR(I+1,I+2))) 

        A3 = -COEF2*(DOT_PROD(I-2,2)*YR(I,I+1) -
     1        DOT_PROD(I-2,2)*YR(I-1,I) + DOT_PROD(I-1,2)*YR(I-2,I-1) +
     1        DOT_PROD(I-1,1)*YR(I-2,I-1) - 2.0*DOT_PROD(I-2,3)*YR(I-1,I) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I-2,1)*YR(I-1,I) -
     1        DOT_PROD(I-2,2)*YR(I-2,I-1)) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I,1)*YR(I-1,I) -
     1        DOT_PROD(I-1,1)*YR(I,I+1) - DOT_PROD(I-1,2)*YR(I,I+1) +
     1        DOT_PROD(I-1,2)*YR(I-1,I))) 

        A4 = -COEF3*(DOT_PROD(I-3,2)*YR(I-2,I-1) -
     1        DOT_PROD(I-2,1)*YR(I-3,I-2) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-2,2)*DOT_PROD(I-3,2) -
     1        DOT_PROD(I-3,3)*DOT_PROD(I-2,1))*(DOT_PROD(I-2,1)*YR(I-1,I) -
     1        DOT_PROD(I-2,2)*YR(I-2,I-1))) 

        FTA_Y(I) = A1 + A2 + A3 + A4 

        A1 = -COEF*(-DOT_PROD(I+1,2)*ZR(I+1,I+2) +
     1        DOT_PROD(I+1,1)*ZR(I+2,I+3) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I+1,2)*DOT_PROD(I,2) -
     1        DOT_PROD(I,3)*DOT_PROD(I+1,1))*(-DOT_PROD(I+1,1)*ZR(I,I+1) +
     1        DOT_PROD(I,2)*ZR(I+1,I+2))) 

        A2 = -COEF1*(-DOT_PROD(I-1,2)*ZR(I+1,I+2) +
     1        DOT_PROD(I,2)*ZR(I,I+1) - DOT_PROD(I,2)*ZR(I-1,I) -
     1        DOT_PROD(I,1)*ZR(I+1,I+2) + 2.0*DOT_PROD(I-1,3)*ZR(I,I+1) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(DOT_PROD(I,1)*ZR(I-1,I) -
     1        DOT_PROD(I-1,1)*ZR(I,I+1) - DOT_PROD(I-1,2)*ZR(I,I+1) +
     1        DOT_PROD(I-1,2)*ZR(I-1,I)) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(-DOT_PROD(I+1,1)*ZR(I,I+1) +
     1        DOT_PROD(I,2)*ZR(I+1,I+2))) 

        A3 = -COEF2*(DOT_PROD(I-2,2)*ZR(I,I+1) -
     1        DOT_PROD(I-2,2)*ZR(I-1,I) + DOT_PROD(I-1,2)*ZR(I-2,I-1) +
     1        DOT_PROD(I-1,1)*ZR(I-2,I-1) - 2.0*DOT_PROD(I-2,3)*ZR(I-1,I) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I-2,1)*ZR(I-1,I) -
     1        DOT_PROD(I-2,2)*ZR(I-2,I-1)) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I,1)*ZR(I-1,I) -
     1        DOT_PROD(I-1,1)*ZR(I,I+1) - DOT_PROD(I-1,2)*ZR(I,I+1) +
     1        DOT_PROD(I-1,2)*ZR(I-1,I))) 
        
        A4 = -COEF3*(DOT_PROD(I-3,2)*ZR(I-2,I-1) -
     1        DOT_PROD(I-2,1)*ZR(I-3,I-2) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-2,2)*DOT_PROD(I-3,2) -
     1        DOT_PROD(I-3,3)*DOT_PROD(I-2,1))*(DOT_PROD(I-2,1)*ZR(I-1,I) -
     1        DOT_PROD(I-2,2)*ZR(I-2,I-1))) 

        FTA_Z(I) = A1 + A2 + A3 + A4 

        ENDDO

C PARTICLE N-2

        I = N-2
        COEF1=(C_PARAM(I)+D_PARAM(I)*(12.0*DCOS(TOR_ANGLE(I))
     1        *DCOS(TOR_ANGLE(I)) - 
     1        3.0))*(1.0D0/DSQRT(X_PROD(I)*X_PROD(I-1)))  

        COEF2=(C_PARAM(I-1)+D_PARAM(I-1)*(12.0*DCOS(TOR_ANGLE(I-1))
     1        *DCOS(TOR_ANGLE(I-1)) - 
     1        3.0))*(1.0D0/DSQRT(X_PROD(I-1)*X_PROD(I-2)))  

        COEF3=(C_PARAM(I-2)+D_PARAM(I-2)*(12.0*DCOS(TOR_ANGLE(I-2))
     1        *DCOS(TOR_ANGLE(I-2)) - 
     1  3.0))*(1.0D0/DSQRT(X_PROD(I-2)*X_PROD(I-3)))  

        A1 = -COEF1*(-DOT_PROD(I-1,2)*XR(I+1,I+2) + 
     1        DOT_PROD(I,2)*XR(I,I+1) - DOT_PROD(I,2)*XR(I-1,I) -
     1        DOT_PROD(I,1)*XR(I+1,I+2) + 2.0*DOT_PROD(I-1,3)*XR(I,I+1) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(DOT_PROD(I,1)*XR(I-1,I) -
     1        DOT_PROD(I-1,1)*XR(I,I+1) - DOT_PROD(I-1,2)*XR(I,I+1) +
     1        DOT_PROD(I-1,2)*XR(I-1,I)) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(-DOT_PROD(I+1,1)*XR(I,I+1) +
     1        DOT_PROD(I,2)*XR(I+1,I+2)))


        A2 = -COEF2*(DOT_PROD(I-2,2)*XR(I,I+1) -
     1        DOT_PROD(I-2,2)*XR(I-1,I) + DOT_PROD(I-1,2)*XR(I-2,I-1) +
     1        DOT_PROD(I-1,1)*XR(I-2,I-1) - 2.0*DOT_PROD(I-2,3)*XR(I-1,I) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I-2,1)*XR(I-1,I) -
     1        DOT_PROD(I-2,2)*XR(I-2,I-1)) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I,1)*XR(I-1,I) -
     1        DOT_PROD(I-1,1)*XR(I,I+1) - DOT_PROD(I-1,2)*XR(I,I+1) +
     1        DOT_PROD(I-1,2)*XR(I-1,I))) 
        
        A3 = -COEF3*(DOT_PROD(I-3,2)*XR(I-2,I-1) -
     1        DOT_PROD(I-2,1)*XR(I-3,I-2) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-2,2)*DOT_PROD(I-3,2) -
     1        DOT_PROD(I-3,3)*DOT_PROD(I-2,1))*(DOT_PROD(I-2,1)*XR(I-1,I) -
     1        DOT_PROD(I-2,2)*XR(I-2,I-1))) 

        FTA_X(I) = A1 + A2 + A3 

        A1 =  -COEF1*(-DOT_PROD(I-1,2)*YR(I+1,I+2) +  
     1        DOT_PROD(I,2)*YR(I,I+1) - DOT_PROD(I,2)*YR(I-1,I) -
     1        DOT_PROD(I,1)*YR(I+1,I+2) + 2.0*DOT_PROD(I-1,3)*YR(I,I+1) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(DOT_PROD(I,1)*YR(I-1,I) -
     1        DOT_PROD(I-1,1)*YR(I,I+1) - DOT_PROD(I-1,2)*YR(I,I+1) +
     1        DOT_PROD(I-1,2)*YR(I-1,I)) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(-DOT_PROD(I+1,1)*YR(I,I+1) +
     1        DOT_PROD(I,2)*YR(I+1,I+2))) 

        A2 =  -COEF2*(DOT_PROD(I-2,2)*YR(I,I+1) -
     1        DOT_PROD(I-2,2)*YR(I-1,I) + DOT_PROD(I-1,2)*YR(I-2,I-1) +
     1        DOT_PROD(I-1,1)*YR(I-2,I-1) - 2.0*DOT_PROD(I-2,3)*YR(I-1,I) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I-2,1)*YR(I-1,I) -
     1        DOT_PROD(I-2,2)*YR(I-2,I-1)) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I,1)*YR(I-1,I) -
     1        DOT_PROD(I-1,1)*YR(I,I+1) - DOT_PROD(I-1,2)*YR(I,I+1) +
     1        DOT_PROD(I-1,2)*YR(I-1,I)))

        A3 = -COEF3*(DOT_PROD(I-3,2)*YR(I-2,I-1) -
     1        DOT_PROD(I-2,1)*YR(I-3,I-2) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-2,2)*DOT_PROD(I-3,2) -
     1        DOT_PROD(I-3,3)*DOT_PROD(I-2,1))*(DOT_PROD(I-2,1)*YR(I-1,I) -
     1        DOT_PROD(I-2,2)*YR(I-2,I-1))) 

        FTA_Y(I) = A1 + A2 + A3 
 
        A1 = -COEF1*(-DOT_PROD(I-1,2)*ZR(I+1,I+2) +  
     1        DOT_PROD(I,2)*ZR(I,I+1) - DOT_PROD(I,2)*ZR(I-1,I) -
     1        DOT_PROD(I,1)*ZR(I+1,I+2) + 2.0*DOT_PROD(I-1,3)*ZR(I,I+1) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(DOT_PROD(I,1)*ZR(I-1,I) -
     1        DOT_PROD(I-1,1)*ZR(I,I+1) - DOT_PROD(I-1,2)*ZR(I,I+1) +
     1        DOT_PROD(I-1,2)*ZR(I-1,I)) -
     1        (1.0D0/X_PROD(I))*(DOT_PROD(I,2)*DOT_PROD(I-1,2) -
     1        DOT_PROD(I-1,3)*DOT_PROD(I,1))*(-DOT_PROD(I+1,1)*ZR(I,I+1) +
     1        DOT_PROD(I,2)*ZR(I+1,I+2))) 

        A2 = -COEF2*(DOT_PROD(I-2,2)*ZR(I,I+1) -
     1        DOT_PROD(I-2,2)*ZR(I-1,I) + DOT_PROD(I-1,2)*ZR(I-2,I-1) +
     1        DOT_PROD(I-1,1)*ZR(I-2,I-1) - 2.0*DOT_PROD(I-2,3)*ZR(I-1,I) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I-2,1)*ZR(I-1,I) -
     1        DOT_PROD(I-2,2)*ZR(I-2,I-1)) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I,1)*ZR(I-1,I) -
     1        DOT_PROD(I-1,1)*ZR(I,I+1) - DOT_PROD(I-1,2)*ZR(I,I+1) +
     1        DOT_PROD(I-1,2)*ZR(I-1,I))) 

        A3 = -COEF3*(DOT_PROD(I-3,2)*ZR(I-2,I-1) -
     1        DOT_PROD(I-2,1)*ZR(I-3,I-2) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-2,2)*DOT_PROD(I-3,2) -
     1        DOT_PROD(I-3,3)*DOT_PROD(I-2,1))*(DOT_PROD(I-2,1)*ZR(I-1,I) -
     1        DOT_PROD(I-2,2)*ZR(I-2,I-1))) 

        FTA_Z(I) = A1 + A2 + A3 

C PARTICLE N-1

        I = N-1
        COEF2=(C_PARAM(I-1)+D_PARAM(I-1)*(12.0*DCOS(TOR_ANGLE(I-1))
     1        *DCOS(TOR_ANGLE(I-1)) - 
     1        3.0))*(1.0D0/DSQRT(X_PROD(I-1)*X_PROD(I-2)))  

        COEF3=(C_PARAM(I-2)+D_PARAM(I-2)*(12.0*DCOS(TOR_ANGLE(I-2))
     1        *DCOS(TOR_ANGLE(I-2)) - 
     1        3.0))*(1.0D0/DSQRT(X_PROD(I-2)*X_PROD(I-3)))  

        A1 = -COEF2*(DOT_PROD(I-2,2)*XR(I,I+1) - 
     1        DOT_PROD(I-2,2)*XR(I-1,I) +
     1        DOT_PROD(I-1,2)*XR(I-2,I-1) +  DOT_PROD(I-1,1)*XR(I-2,I-1) -
     1        2.0*DOT_PROD(I-2,3)*XR(I-1,I) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I-2,1)*XR(I-1,I) -
     1        DOT_PROD(I-2,2)*XR(I-2,I-1)) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I,1)*XR(I-1,I) -
     1        DOT_PROD(I-1,1)*XR(I,I+1) - DOT_PROD(I-1,2)*XR(I,I+1) +
     1        DOT_PROD(I-1,2)*XR(I-1,I))) 

        A2 = -COEF3*(DOT_PROD(I-3,2)*XR(I-2,I-1) - 
     1        DOT_PROD(I-2,1)*XR(I-3,I-2) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-2,2)*DOT_PROD(I-3,2) -
     1        DOT_PROD(I-3,3)*DOT_PROD(I-2,1))*(DOT_PROD(I-2,1)*XR(I-1,I) -
     1        DOT_PROD(I-2,2)*XR(I-2,I-1))) 

        FTA_X(I) = A1 + A2  

        A1 = -COEF2*(DOT_PROD(I-2,2)*YR(I,I+1) - 
     1        DOT_PROD(I-2,2)*YR(I-1,I) +
     1        DOT_PROD(I-1,2)*YR(I-2,I-1) +  DOT_PROD(I-1,1)*YR(I-2,I-1) -
     1        2.0*DOT_PROD(I-2,3)*YR(I-1,I) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I-2,1)*YR(I-1,I) -
     1        DOT_PROD(I-2,2)*YR(I-2,I-1)) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I,1)*YR(I-1,I) -
     1  DOT_PROD(I-1,1)*YR(I,I+1) - DOT_PROD(I-1,2)*YR(I,I+1) +
     1        DOT_PROD(I-1,2)*YR(I-1,I))) 

        A2 = -COEF3*(DOT_PROD(I-3,2)*YR(I-2,I-1) - 
     1        DOT_PROD(I-2,1)*YR(I-3,I-2) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-2,2)*DOT_PROD(I-3,2) -
     1        DOT_PROD(I-3,3)*DOT_PROD(I-2,1))*(DOT_PROD(I-2,1)*YR(I-1,I) -
     1        DOT_PROD(I-2,2)*YR(I-2,I-1))) 

        FTA_Y(I) = A1 + A2  

        A1 = -COEF2*(DOT_PROD(I-2,2)*ZR(I,I+1) - 
     1        DOT_PROD(I-2,2)*ZR(I-1,I) +
     1        DOT_PROD(I-1,2)*ZR(I-2,I-1) +  DOT_PROD(I-1,1)*ZR(I-2,I-1) -
     1        2.0*DOT_PROD(I-2,3)*ZR(I-1,I) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I-2,1)*ZR(I-1,I) -
     1        DOT_PROD(I-2,2)*ZR(I-2,I-1)) -
     1        (1.0D0/X_PROD(I-1))*(DOT_PROD(I-1,2)*DOT_PROD(I-2,2) -
     1        DOT_PROD(I-2,3)*DOT_PROD(I-1,1))*(DOT_PROD(I,1)*ZR(I-1,I) -
     1        DOT_PROD(I-1,1)*ZR(I,I+1) - DOT_PROD(I-1,2)*ZR(I,I+1) +
     1        DOT_PROD(I-1,2)*ZR(I-1,I))) 

        A2 = -COEF3*(DOT_PROD(I-3,2)*ZR(I-2,I-1) - 
     1        DOT_PROD(I-2,1)*ZR(I-3,I-2) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-2,2)*DOT_PROD(I-3,2) -
     1        DOT_PROD(I-3,3)*DOT_PROD(I-2,1))*(DOT_PROD(I-2,1)*ZR(I-1,I) -
     1        DOT_PROD(I-2,2)*ZR(I-2,I-1))) 

        FTA_Z(I) = A1 + A2 
 
C PARTICLE N

        I = N
        COEF3=(C_PARAM(I-2)+D_PARAM(I-2)*(12.0*DCOS(TOR_ANGLE(I-2))
     1        *DCOS(TOR_ANGLE(I-2)) - 
     1        3.0))*(1.0D0/DSQRT(X_PROD(I-2)*X_PROD(I-3)))  

        FTA_X(I) = -COEF3*(DOT_PROD(I-3,2)*XR(I-2,I-1) 
     1        - DOT_PROD(I-2,1)*XR(I-3,I-2) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-2,2)*DOT_PROD(I-3,2) -
     1        DOT_PROD(I-3,3)*DOT_PROD(I-2,1))*(DOT_PROD(I-2,1)*XR(I-1,I) -
     1        DOT_PROD(I-2,2)*XR(I-2,I-1))) 

        FTA_Y(I) = -COEF3*(DOT_PROD(I-3,2)*YR(I-2,I-1) - 
     1        DOT_PROD(I-2,1)*YR(I-3,I-2) 
     1        - (1.0D0/X_PROD(I-2))*(DOT_PROD(I-2,2)*DOT_PROD(I-3,2) 
     1        - DOT_PROD(I-3,3)*DOT_PROD(I-2,1))*(DOT_PROD(I-2,1)*YR(I-1,I) -
     1        DOT_PROD(I-2,2)*YR(I-2,I-1))) 

        FTA_Z(I) = -COEF3*(DOT_PROD(I-3,2)*ZR(I-2,I-1) - 
     1        DOT_PROD(I-2,1)*ZR(I-3,I-2) -
     1        (1.0D0/X_PROD(I-2))*(DOT_PROD(I-2,2)*DOT_PROD(I-3,2) -
     1        DOT_PROD(I-3,3)*DOT_PROD(I-2,1))*(DOT_PROD(I-2,1)*ZR(I-1,I) -
     1        DOT_PROD(I-2,2)*ZR(I-2,I-1))) 

C TOTAL UP THE GRADIENTS

        DO I = 1, N
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

        RETURN
        END
C }}}
C> CALCULATE THE SECOND DERIVATIVE MATRIX (TWO-SIDED NUMERICAL APPROACH)

        SUBROUTINE CALC_DYN(QO,N,A_PARAM,B_PARAM,C_PARAM,D_PARAM,X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD, 
     1                      BOND_ANGLE,TOR_ANGLE,RADII,NTYPE)
C {{{
C DECLARATIONS {{{
        USE MODHESS
        IMPLICIT NONE
        DOUBLE PRECISION RMASS, EPSILON, SIGMA, DELTA, THETA_0, RK_R, RK_THETA
        PARAMETER (RMASS = 40.0, EPSILON = 0.0100570)
        PARAMETER (SIGMA=3.4 ,DELTA=1.0D-4, THETA_0 = 1.8326)
        PARAMETER (RK_R = 20.0*0.0100570, RK_THETA = 20.0*0.0100570)
        INTEGER NTYPE(46), N, I, J
        DOUBLE PRECISION QO(3*N), FQ1(3*N), 
     1  FQ2(3*N)
        DOUBLE PRECISION A_PARAM(N,N), B_PARAM(N,N), D_PARAM(N),C_PARAM(N),
     1                  X(N), Y(N),Z(N),XR(N,N),YR(N,N),ZR(N,N),
     2                  DOT_PROD(N,3), X_PROD(N), BOND_ANGLE(N), TOR_ANGLE(N), 
     3                  RADII(N,N)
C }}}
C FILL IN THE HESSIAN MATRIX
C {{{

        DO J = 1, 3*N

        QO(J) = QO(J) + DELTA
        CALL CALC_INT_COORDS(QO,N,A_PARAM,B_PARAM,C_PARAM,D_PARAM,X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD,
     1                       BOND_ANGLE,TOR_ANGLE,RADII,NTYPE)
        CALL CALC_GRADIENT(QO,FQ2,N,A_PARAM,B_PARAM,C_PARAM,D_PARAM,X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD,
     1                     BOND_ANGLE,TOR_ANGLE,RADII,NTYPE)
        QO(J) = QO(J) - 2.0*DELTA
        CALL CALC_INT_COORDS(QO,N,A_PARAM,B_PARAM,C_PARAM,D_PARAM,X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD,
     1                       BOND_ANGLE,TOR_ANGLE,RADII,NTYPE)
        CALL CALC_GRADIENT(QO,FQ1,N,A_PARAM,B_PARAM,C_PARAM,D_PARAM,X,Y,Z,XR,YR,ZR,DOT_PROD,X_PROD,
     1                     BOND_ANGLE,TOR_ANGLE,RADII,NTYPE)
        QO(J) = QO(J) + DELTA

        DO I = J, 3*N

        HESS(I,J) = (FQ2(I) -  FQ1(I))/(2.0*DELTA)
        HESS(J,I) = HESS(I,J)

        ENDDO
        ENDDO

C }}}
        RETURN
        END
C }}}
C> FILL THE PARAMETER ARRAYS

        SUBROUTINE PARAM_ARRAY(A_PARAM,B_PARAM,C_PARAM,D_PARAM,N)
C {{{
C DECLARATIONS {{{
        IMPLICIT NONE
        INTEGER NTYPE(46), N, ICOUNT, J, I
        DOUBLE PRECISION A_PARAM(N,N), B_PARAM(N,N), EPSILON
        DOUBLE PRECISION C_PARAM(N), D_PARAM(N)
        PARAMETER (EPSILON = 0.0100570D0)
C }}}
C FIRSTLY, SPECIFY AMINO ACID TYPES BY FILLING IN ARRAY NTYPE(:)
C 1 -> HYDROPHOBIC (B)
C 2 -> HYDROPHILIC (L)
C 3 -> NEUTRAL (N)
C {{{
        NTYPE(1) = 1
        NTYPE(2) = 1
        NTYPE(3) = 1
        NTYPE(4) = 1
        NTYPE(5) = 1
        NTYPE(6) = 1
        NTYPE(7) = 1
        NTYPE(8) = 1
        NTYPE(9) = 1
        NTYPE(10) = 3
        NTYPE(11) = 3
        NTYPE(12) = 3
        NTYPE(13) = 2
        NTYPE(14) = 1
        NTYPE(15) = 2
        NTYPE(16) = 1
        NTYPE(17) = 2
        NTYPE(18) = 1
        NTYPE(19) = 2
        NTYPE(20) = 1
        NTYPE(21) = 3
        NTYPE(22) = 3
        NTYPE(23) = 3
        NTYPE(24) = 1
        NTYPE(25) = 1
        NTYPE(26) = 1
        NTYPE(27) = 1
        NTYPE(28) = 1
        NTYPE(29) = 1
        NTYPE(30) = 1
        NTYPE(31) = 1
        NTYPE(32) = 1
        NTYPE(33) = 3
        NTYPE(34) = 3
        NTYPE(35) = 3
        NTYPE(36) = 2
        NTYPE(37) = 1
        NTYPE(38) = 2
        NTYPE(39) = 1
        NTYPE(40) = 2
        NTYPE(41) = 1
        NTYPE(42) = 2
        NTYPE(43) = 1
        NTYPE(44) = 2
        NTYPE(45) = 1
        NTYPE(46) = 2
C }}}
C SPECIFY PARAMETERS FOR THE DIHEDRAL ANGLE POTENTIAL.
C       THESE ARE STORED IN ARRAYS 
C       C_PARAM(:), D_PARAM(:) 
C {{{

        DO I = 1, N-3
        ICOUNT = 0

        DO J = 0,3
        IF(NTYPE(I+J) .EQ. 3)THEN
        ICOUNT = ICOUNT + 1
        ENDIF
        ENDDO

        IF(ICOUNT .GE. 2)THEN
        C_PARAM(I+1) = 0.0
        D_PARAM(I+1) = 0.2*EPSILON
        ELSE
        C_PARAM(I+1) = 1.2*EPSILON
        D_PARAM(I+1) = 1.2*EPSILON
        ENDIF

        ICOUNT = 0

        ENDDO
!}}}
C SPECIFY PARAMETERS FOR THE L-J INTERACTION BETWEEN NON-BONDED PARTICLES.
C       THESE ARE STORED IN ARRAYS 
C       A_PARAM(:,:), B_PARAM(:,:) 
C {{{

        DO I = 1, N-1
        DO J = I+1, N

        IF (NTYPE(I) .EQ. 3 .OR. NTYPE(J) .EQ. 3)THEN
        A_PARAM(I,J) = 1.0*EPSILON 
        B_PARAM(I,J) = 0.0 
        A_PARAM(J,I) = 1.0*EPSILON 
        B_PARAM(J,I) = 0.0

        ELSEIF (NTYPE(I) .EQ. 1 .AND. NTYPE(J) .EQ. 1)THEN
        A_PARAM(I,J) =  EPSILON
        B_PARAM(I,J) = -EPSILON 
        A_PARAM(J,I) =  EPSILON
        B_PARAM(J,I) = -EPSILON
        
        ELSE

        A_PARAM(I,J) = EPSILON*2.0/3.0 
        B_PARAM(I,J) = EPSILON*2.0/3.0 
        A_PARAM(J,I) = EPSILON*2.0/3.0 
        B_PARAM(J,I) = EPSILON*2.0/3.0 

        ENDIF

        ENDDO
        ENDDO
!}}}
        RETURN
        END
C }}}
