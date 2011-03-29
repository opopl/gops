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
C
C********************************************************************
C
C SUBROUTINE SCDERIV CALCULATES FIRST AND SECOND DERIVATIVE 
C ANALYTICALLY:
C
C********************************************************************
C
      SUBROUTINE SCDIFF(N,X,V,EPS,C,SIG,NN,MM,PSC)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6, J7, NN, MM, I, J
      
      DOUBLE PRECISION X(3*N), 
     1                 V(3*N), RHO(N), RM(N,N),
     2                 RN(N,N),RN2(N,N),
     3                 RM2(N,N), RN4(N,N),EPS,C,
     4                 RM4(N,N), SIG,PSC,DIST,POTA,POTB,TEMP,TEMP1,TEMP2,TEMP3
C
C
C  STORE DISTANCE MATRICES.
C
      DO 20 J1=1,N
         RM(J1,J1)=0.0D0
         RN(J1,J1)=0.0D0
         RM2(J1,J1)=0.0D0
         RM4(J1,J1)=0.0D0
         RN2(J1,J1)=0.0D0
         RN4(J1,J1)=0.0D0
         DO 10 J2=1,J1-1
            DIST=(X(3*(J1-1)+1)-X(3*(J2-1)+1))**2
     1          +(X(3*(J1-1)+2)-X(3*(J2-1)+2))**2
     2          +(X(3*(J1-1)+3)-X(3*(J2-1)+3))**2
            RM(J2,J1)=DIST**(-FLOAT(MM)/2.0D0)
            RN(J2,J1)=DIST**(-FLOAT(NN)/2.0D0)
            RM2(J2,J1)=DIST**(-(FLOAT(MM)+2.0D0)/2.0D0)
            RM4(J2,J1)=DIST**(-(FLOAT(MM)+4.0D0)/2.0D0)
            RN2(J2,J1)=DIST**(-(FLOAT(NN)+2.0D0)/2.0D0)
            RN4(J2,J1)=DIST**(-(FLOAT(NN)+4.0D0)/2.0D0)
            RM(J1,J2)=RM(J2,J1)
            RN(J1,J2)=RN(J2,J1)
            RM2(J1,J2)=RM2(J2,J1)
            RM4(J1,J2)=RM4(J2,J1)
            RN2(J1,J2)=RN2(J2,J1)
            RN4(J1,J2)=RN4(J2,J1)
10       CONTINUE
20    CONTINUE
C
C STORE DENSITY MATRIX: IN THE CASE OF THE PERFECT FCC LATTICE,
C THE INFINITELY EXTENDED CRYSTAL IMPLIES THAT EVERY RHO(J) IS
C EQUAL TO RHO(1).
C
      DO 11 I=1,N
         RHO(I)=0.0D00
         DO 122 J=1,N
            RHO(I)=RHO(I) + SIG**MM*RM(I,J)
122      CONTINUE
11    CONTINUE
C
C FIRST CALCULATE THE POTENTIAL ENERGY:
C
      POTA=0.0D0
      POTB=0.0D0
      DO 13 I=1,N
         DO 14 J=1,N
            POTA=POTA + 0.50D00*EPS*SIG**NN*RN(I,J)
14       CONTINUE
         POTB=POTB + EPS*DSQRT(RHO(I))*C
13    CONTINUE
      PSC=POTA - POTB
      PRINT*, ' SUTTON-CHEN N=', NN,' M=', MM
      PRINT*,'TWO-BODY CONTRIBUTION=',POTA,' EV'
      PRINT*,'MANY-BODY CONTRIBUTION=',-POTB,' EV'
C     PRINT*, 'TOTAL ENERGY FOR LAST STEP=', PSC,' EV'
C
C NOW CALCULATE THE GRADIENT ANALYTICALLY.
C
      DO 50 J1=1,N
         DO 40 J2=1,3
            J3=3*(J1-1)+J2
            V(J3)=0.0D0
            DO 30 J4=1,N
               V(J3)=V(J3)+( -NN*EPS*SIG**NN*RN2(J1,J4)
     1                    +FLOAT(MM)/2*EPS*C*SIG**MM*RM2(J1,J4)*
     2                     (1/DSQRT(RHO(J1)) + 1/DSQRT(RHO(J4))))
     3              *(X(J3)-X(3*(J4-1)+J2))
30          CONTINUE
C           PRINT*,'V(',J3,')=',V(J3)
40       CONTINUE
50    CONTINUE
C
C  NOW DO THE HESSIAN. FIRST ARE THE ENTIRELY DIAGONAL TERMS.
C
      DO 80 J1=1,N
         DO 70 J2=1,3
            J3=3*(J1-1)+J2
            HESS(J3,J3)=0.0D0
            TEMP=0.0D0
            DO 60 J4=1,N
               HESS(J3,J3)=HESS(J3,J3)
     1             - NN*EPS*(
     2                     SIG**NN*RN2(J4,J1)
     3                   - SIG**NN*RN4(J4,J1)*(NN+2)
     4                   * (X(J3)-X(3*(J4-1)+J2))**2 )
     5             + FLOAT(MM)*EPS*C/2*(
     6                   ( SIG**MM*RM2(J4,J1) 
     7                   - SIG**MM*RM4(J4,J1)*(MM+2)
     8                   * (X(J3)-X(3*(J4-1)+J2))**2 )
     9                   *(1/DSQRT(RHO(J4)) + 1/DSQRT(RHO(J1)))
     A                   + 0.5*SIG**(2*MM)*MM*RM2(J4,J1)**2
     B                   *(X(J3)-X(3*(J4-1)+J2))**2
     C                   *RHO(J4)**(-1.5) )
               TEMP=TEMP+(X(J3)-X(3*(J4-1)+J2))*RM2(J4,J1)
60          CONTINUE
            TEMP=MM*MM*EPS*C*SIG**(2*MM)*TEMP*TEMP/(4.0*RHO(J1)**(1.5)) 
            HESS(J3,J3)=HESS(J3,J3)+TEMP
C           PRINT*,'HESS(',J3,',',J3,')=',HESS(J3,J3)
70       CONTINUE
80    CONTINUE
C
C  NEXT ARE THE TERMS WHERE X_I AND X_J ARE ON THE SAME ATOM
C  BUT ARE DIFFERENT, E.G. Y AND Z.
C
      DO 120 J1=1,N
         DO 110 J2=1,3
            J3=3*(J1-1)+J2
            DO 100 J4=1,J2-1
               HESS(J3,3*(J1-1)+J4)=0.0D0
               TEMP1=0.0D0
               TEMP2=0.0D0
               DO 90 J5=1,N
                  HESS(J3,3*(J1-1)+J4)=HESS(J3,3*(J1-1)+J4)
     1           + NN*EPS*(SIG**NN*RN4(J1,J5)*(NN+2)
     2          *(X(J3)-X(3*(J5-1)+J2))*(X(3*(J1-1)+J4)-X(3*(J5-1)+J4)))
     3          +  FLOAT(MM)*EPS*C/2*(
     4            -(1/DSQRT(RHO(J5)) + 1/DSQRT(RHO(J1)))
     5            * SIG**MM*RM4(J5,J1)*(MM+2)
     6           *(X(J3)-X(3*(J5-1)+J2))*(X(3*(J1-1)+J4)-X(3*(J5-1)+J4))
     7           + 0.5*SIG**(2*MM)*MM*RM2(J5,J1)**2
     8           *(X(J3)-X(3*(J5-1)+J2))*(X(3*(J1-1)+J4)-X(3*(J5-1)+J4))
     9           *RHO(J5)**(-1.5)  )
                  TEMP1=TEMP1+(X(J3)-X(3*(J5-1)+J2))*RM2(J5,J1)
                  TEMP2=TEMP2+(X(3*(J1-1)+J4)-X(3*(J5-1)+J4))*RM2(J5,J1)
90             CONTINUE
               TEMP=MM*MM*EPS*C*SIG**(2*MM)*TEMP1*TEMP2/
     1              (4.0*RHO(J1)**(1.5))  
               HESS(J3,3*(J1-1)+J4)=HESS(J3,3*(J1-1)+J4)+TEMP
               HESS(3*(J1-1)+J4,J3)=HESS(J3,3*(J1-1)+J4)
C              PRINT*,'HESS(',J3,',',3*(J1-1)+J4,')=',HESS(J3,3*(J1-1)+J4)
100         CONTINUE
110      CONTINUE
120   CONTINUE
C
C  CASE III, DIFFERENT ATOMS, SAME CARTESIAN COORDINATE.
C
      DO 150 J1=1,N
         DO 140 J2=1,3
            J3=3*(J1-1)+J2
            DO 130 J4=J1+1,N
               HESS(J3,3*(J4-1)+J2)=
     1               NN*EPS*(
     2                     SIG**NN*RN2(J4,J1)
     3                   - SIG**NN*RN4(J4,J1)*(NN+2)
     4                   * (X(J3)-X(3*(J4-1)+J2))**2 )
     5             - FLOAT(MM)*EPS*C/2*(
     6                   ( SIG**MM*RM2(J4,J1) 
     7                   - SIG**MM*RM4(J4,J1)*(MM+2)
     8                   * (X(J3)-X(3*(J4-1)+J2))**2 )
     9                   *(1/DSQRT(RHO(J4)) + 1/DSQRT(RHO(J1))))
               TEMP1=0.0D0
               TEMP2=0.0D0
               TEMP3=0.0D0
               DO 125 J5=1,N
                  TEMP1=TEMP1+(X(3*(J4-1)+J2)-X(3*(J5-1)+J2))*RM2(J5,J4) 
                  TEMP2=TEMP2+(X(J3)-X(3*(J5-1)+J2))*RM2(J5,J1) 
                  TEMP3=TEMP3+(X(3*(J4-1)+J2)-X(3*(J5-1)+J2))*RM2(J5,J4) 
     1                       *(X(J3)-X(3*(J5-1)+J2))*RM2(J5,J1)*   
     2                        RHO(J5)**(-1.5)
125            CONTINUE
              TEMP1=MM*MM*EPS*C*SIG**(2*MM)*TEMP1*(X(J3)-X(3*(J4-1)+J2))  
     1              *RM2(J4,J1)/(4.0*RHO(J4)**(1.5))  
              TEMP2=MM*MM*EPS*C*SIG**(2*MM)*TEMP2*(X(J3)-X(3*(J4-1)+J2)) 
     1              *RM2(J4,J1)/(4.0*RHO(J1)**(1.5))  
               TEMP3=MM*MM*EPS*C*SIG**(2*MM)*TEMP3/4.0
               HESS(J3,3*(J4-1)+J2)=HESS(J3,3*(J4-1)+J2)+TEMP1-TEMP2+TEMP3
               HESS(3*(J4-1)+J2,J3)=HESS(J3,3*(J4-1)+J2)
C              PRINT*,'HESS(',J3,',',3*(J4-1)+J2,')=',HESS(J3,3*(J4-1)+J2)
130         CONTINUE
140      CONTINUE
150   CONTINUE
C
C  CASE IV: DIFFERENT ATOMS AND DIFFERENT CARTESIAN COORDINATES.
C
      DO 180 J1=1,N
         DO 170 J2=1,3
            J3=3*(J1-1)+J2
            DO 160 J4=J1+1,N
               DO 155 J5=1,J2-1
                  J6=3*(J4-1)+J5
                  HESS(J3,J6)=
     1          -  NN*EPS*SIG**NN*RN4(J4,J1)*(NN+2)
     2           *(X(J3)-X(3*(J4-1)+J2))*(X(3*(J1-1)+J5)-X(J6))
     3          +  FLOAT(MM)*EPS*C/2*
     4             (1/DSQRT(RHO(J4)) + 1/DSQRT(RHO(J1)))
     5            * SIG**MM*RM4(J4,J1)*(MM+2)
     6           *(X(J3)-X(3*(J4-1)+J2))*(X(3*(J1-1)+J5)-X(J6))
                  TEMP1=0.0D0
                  TEMP2=0.0D0
                  TEMP3=0.0D0
                  DO 151 J7=1,N
                     TEMP1=TEMP1+(X(J6)-X(3*(J7-1)+J5))*RM2(J7,J4)
                     TEMP2=TEMP2+(X(J3)-X(3*(J7-1)+J2))*RM2(J7,J1)
                     TEMP3=TEMP3+(X(J3)-X(3*(J7-1)+J2))*RM2(J7,J1)*
     1                           (X(J6)-X(3*(J7-1)+J5))*RM2(J7,J4)*   
     2                            RHO(J7)**(-1.5)
151               CONTINUE
                  TEMP1=MM*MM*EPS*C*SIG**(2*MM)*TEMP1*
     1            (X(J3)-X(3*(J4-1)+J2))*RM2(J4,J1)/(4.0*RHO(J4)**(1.5)) 
                  TEMP2=MM*MM*EPS*C*SIG**(2*MM)*TEMP2*
     1            (X(3*(J1-1)+J5)-X(J6))*RM2(J4,J1)/(4.0*RHO(J1)**(1.5)) 
                  TEMP3=MM*MM*EPS*C*SIG**(2*MM)*TEMP3/4.0
                  HESS(J3,J6)=HESS(J3,J6)+TEMP1-TEMP2+TEMP3
                  HESS(J6,J3)=HESS(J3,J6)
C              PRINT*,'HESS(',J3,',',J6,')=',HESS(J3,J6)
155            CONTINUE
               DO 156 J5=J2+1,3
                  J6=3*(J4-1)+J5
                  HESS(J3,J6)=
     1          -  NN*EPS*(SIG**NN*RN4(J1,J4)*(NN+2)
     2           *(X(J3)-X(3*(J4-1)+J2))*(X(3*(J1-1)+J5)-X(J6)))
     3          -  FLOAT(MM)*EPS*C/2*(
     4            -(1/DSQRT(RHO(J4)) + 1/DSQRT(RHO(J1)))
     5            * SIG**MM*RM4(J4,J1)*(MM+2)
     6           *(X(J3)-X(3*(J4-1)+J2))*(X(3*(J1-1)+J5)-X(J6)))
                  TEMP1=0.0D0
                  TEMP2=0.0D0
                  TEMP3=0.0D0
                  DO 152 J7=1,N
                     TEMP1=TEMP1+(X(J6)-X(3*(J7-1)+J5))*RM2(J7,J4)
                     TEMP2=TEMP2+(X(J3)-X(3*(J7-1)+J2))*RM2(J7,J1)
                     TEMP3=TEMP3+(X(J3)-X(3*(J7-1)+J2))*RM2(J7,J1)*
     1                           (X(J6)-X(3*(J7-1)+J5))*RM2(J7,J4)*  
     1                            RHO(J7)**(-1.5)
152               CONTINUE
                  TEMP1=MM*MM*EPS*C*SIG**(2*MM)*TEMP1*
     1            (X(J3)-X(3*(J4-1)+J2))*RM2(J4,J1)/(4.0*RHO(J4)**(1.5)) 
                  TEMP2=MM*MM*EPS*C*SIG**(2*MM)*TEMP2*
     1            (X(3*(J1-1)+J5)-X(J6))*RM2(J4,J1)/(4.0*RHO(J1)**(1.5)) 
                  TEMP3=MM*MM*EPS*C*SIG**(2*MM)*TEMP3/4.0
                  HESS(J3,J6)=HESS(J3,J6)+TEMP1-TEMP2+TEMP3
                  HESS(J6,J3)=HESS(J3,J6)
156            CONTINUE
160         CONTINUE
170      CONTINUE
180   CONTINUE
      RETURN
      END
