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
C       POTENTIAL ENERGY OF A (C60)N SYSTEM
C       PACHECO-PRATES-RAMALHO POTENTIAL (PRL 1997)
C
C      X,Y,Z: VECTORS

      SUBROUTINE PRC60(NATOMS,LCOORDS,V,EPPR,GTEST,SECT)
      USE MODHESS
      IMPLICIT NONE
        DOUBLE PRECISION CAT
C       PARAMETER(CAT=4752000.D0)
        PARAMETER(CAT=0.0D0)
        INTEGER J1, J2, I, J, NATOMS
      DOUBLE PRECISION X(NATOMS),Y(NATOMS),Z(NATOMS),LCOORDS(3*NATOMS),V3B,V2B,EPPR,V(3*NATOMS)
      DOUBLE PRECISION MIJ,FIJ,WIJ,XIJ,YIJ,ZIJ,RIJ,RIJ2,VIJ,PMORSE,WTOT,RIJ23,
     1                   FX(NATOMS),FY(NATOMS),FZ(NATOMS),ERMI,DWIJ,DFIJ,DMIJ
        DOUBLE PRECISION DMU, DELTA
      PARAMETER(DMU=10.05D0)
      PARAMETER(DELTA=1.04D0)
        DOUBLE PRECISION DM0, TAU, D0
      PARAMETER(DM0=0.3D0)
      PARAMETER(TAU=9.75D0)
      PARAMETER(D0=10.3D0)
        DOUBLE PRECISION C6, C8, C10, C12
      PARAMETER(C6=75600.D0)
      PARAMETER(C8=9122400.D0)
      PARAMETER(C10=2.09D8)
      PARAMETER(C12=7.78D10)
        LOGICAL GTEST,SECT

        DO J1=1,NATOMS
           J2=3*(J1-1)
           X(J1)=LCOORDS(J2+1)
           Y(J1)=LCOORDS(J2+2)
           Z(J1)=LCOORDS(J2+3)
        ENDDO

      V3B=0.D0
      V2B=0.D0

        IF (.NOT.GTEST) THEN
      DO I=1,NATOMS
         DO J=I+1,NATOMS
C
C       2-BODY INTERACTION
C
            XIJ=X(J)-X(I)
            YIJ=Y(J)-Y(I)
            ZIJ=Z(J)-Z(I)
            RIJ2=XIJ*XIJ+YIJ*YIJ+ZIJ*ZIJ
            RIJ=SQRT(RIJ2)
            
            FIJ=1.0D0/(1.D0+EXP((RIJ-DMU)/DELTA))
              PMORSE=EXP(TAU*(1.D0-RIJ/D0))
            MIJ=DM0*PMORSE*(PMORSE-2.D0)
            WIJ=-(C6+(C8+(C10+C12/RIJ2)/RIJ2)/RIJ2)/RIJ2**3
            
            VIJ=FIJ*MIJ+(1.D0-FIJ)*WIJ
            
            V2B=V2B+VIJ
         ENDDO
      ENDDO

        ELSE

      DO I=1,NATOMS
           FX(I)=0.D0
           FY(I)=0.D0
           FZ(I)=0.D0
        ENDDO
      DO I=1,NATOMS
         DO J=I+1,NATOMS
C
C       2-BODY INTERACTION
C
            XIJ=X(J)-X(I)
            YIJ=Y(J)-Y(I)
            ZIJ=Z(J)-Z(I)
            RIJ2=XIJ*XIJ+YIJ*YIJ+ZIJ*ZIJ
            RIJ=SQRT(RIJ2)
            
              ERMI=EXP((RIJ-DMU)/DELTA)
              FIJ=1.0D0/(1.D0+ERMI)
              DFIJ=-ERMI/(DELTA*(1.D0+ERMI)**2)

              PMORSE=EXP(TAU*(1.D0-RIJ/D0))
            MIJ=DM0*PMORSE*(PMORSE-2.D0)
              DMIJ=(2.D0*TAU*DM0*PMORSE*(1.D0-PMORSE))/D0

              RIJ23=RIJ2**3
            WIJ=-(C6+(C8+(C10+C12/RIJ2)/RIJ2)/RIJ2)/RIJ23
              DWIJ=(6*C6+(8*C8+(10*C10+12*C12/RIJ2)/RIJ2)/RIJ2)/(RIJ*RIJ23)
            
            VIJ=FIJ*MIJ+(1.D0-FIJ)*WIJ
            
            V2B=V2B+VIJ

              WTOT=MIJ*DFIJ+FIJ*DMIJ+(1.D0-FIJ)*DWIJ-DFIJ*WIJ
              WTOT=-WTOT

              FX(I)=FX(I)+(X(I)-X(J))*WTOT/RIJ
              FY(I)=FY(I)+(Y(I)-Y(J))*WTOT/RIJ
              FZ(I)=FZ(I)+(Z(I)-Z(J))*WTOT/RIJ

              FX(J)=FX(J)+(X(J)-X(I))*WTOT/RIJ
              FY(J)=FY(J)+(Y(J)-Y(I))*WTOT/RIJ
              FZ(J)=FZ(J)+(Z(J)-Z(I))*WTOT/RIJ
       
         ENDDO
      ENDDO
        ENDIF
      
      EPPR=V2B

        IF (GTEST) THEN
           DO J1=1,NATOMS
              J2=3*(J1-1)
              V(J2+1)=-FX(J1)
              V(J2+2)=-FY(J1)
              V(J2+3)=-FZ(J1)
C             WRITE(*,'(A,I3,3F20.10)') 'J1,FX,FY,FZ=',J1,FX(J1),FY(J1),FZ(J1)
           ENDDO
        ENDIF

        IF (SECT) CALL HESSIAN(X,Y,Z,NATOMS)
      
      RETURN
      END

      SUBROUTINE HESSIAN(X,Y,Z,NSIZE)
      USE MODHESS
      IMPLICIT NONE
        INTEGER NSIZE,I,J
      DOUBLE PRECISION X(NSIZE),Y(NSIZE),Z(NSIZE),WYZ,WXZ,WXY,WZZ,WYY,WXX,D2V,D2V2,DV2,DVDW,VDW
      DOUBLE PRECISION MIJ,DMIJ,FIJ,DFIJ,WIJ,DWIJ,XMORSE,XX,YY,ZZ,XY,XZ,YZ,XIJ,YIJ,ZIJ,RIJ2,RIJ,FERMI,DFERMI,DMORSE
      DOUBLE PRECISION D2MIJ,D2FIJ,D2WIJ,D2MORSE,D2FERMI,D2VDW

      DO I=1,NSIZE
         XX=0.0D0
         YY=0.0D0
         ZZ=0.0D0
         XY=0.0D0
         XZ=0.0D0
         YZ=0.0D0
         DO J=1,NSIZE
            IF (I.NE.J) THEN
         XIJ=X(J)-X(I)
         YIJ=Y(J)-Y(I)
         ZIJ=Z(J)-Z(I)
         RIJ2=XIJ*XIJ+YIJ*YIJ+ZIJ*ZIJ
         RIJ=DSQRT(RIJ2)
         
         FIJ=FERMI(RIJ)
         DFIJ=DFERMI(RIJ)
         D2FIJ=D2FERMI(RIJ)
         MIJ=XMORSE(RIJ)
         DMIJ=DMORSE(RIJ)
         D2MIJ=D2MORSE(RIJ)
         WIJ=VDW(RIJ)
         DWIJ=DVDW(RIJ)
         D2WIJ=D2VDW(RIJ)
         
         DV2=MIJ*DFIJ+FIJ*DMIJ+(1.D0-FIJ)*DWIJ-DFIJ*WIJ
         D2V2=D2FIJ*MIJ+2.D0*DFIJ*DMIJ+D2MIJ*FIJ
         D2V2=D2V2+(1.D0-FIJ)*D2WIJ-D2FIJ*WIJ-2.D0*DFIJ*DWIJ

         D2V=DV2/RIJ-D2V2

         WXX=-DV2/RIJ+XIJ**2*D2V/RIJ2
         WYY=-DV2/RIJ+YIJ**2*D2V/RIJ2
         WZZ=-DV2/RIJ+ZIJ**2*D2V/RIJ2
         WXY=XIJ*YIJ*D2V/RIJ2
         WXZ=XIJ*ZIJ*D2V/RIJ2
         WYZ=YIJ*ZIJ*D2V/RIJ2

         HESS(3*I-2,3*J-2)=WXX
         HESS(3*I-1,3*J-1)=WYY
         HESS(3*I,3*J)    =WZZ
         HESS(3*I-2,3*J-1)=WXY
         HESS(3*I-2,3*J)=WXZ
         HESS(3*I-1,3*J)=WYZ
         HESS(3*I-1,3*J-2)=WXY
         HESS(3*I,3*J-2)=WXZ
         HESS(3*I,3*J-1)=WYZ

         XX=XX+WXX
         YY=YY+WYY
         ZZ=ZZ+WZZ
         XY=XY+WXY
         XZ=XZ+WXZ
         YZ=YZ+WYZ

            ENDIF

         ENDDO

         HESS(3*I-2,3*I-2)=-XX
         HESS(3*I-1,3*I-1)=-YY
         HESS(3*I,3*I)=-ZZ
         HESS(3*I-2,3*I-1)=-XY
         HESS(3*I-2,3*I)=-XZ
         HESS(3*I-1,3*I)=-YZ
         HESS(3*I-1,3*I-2)=-XY
         HESS(3*I,3*I-2)=-XZ
         HESS(3*I,3*I-1)=-YZ
      
      ENDDO

      RETURN
      END

C__________________________________________________________________________

      FUNCTION FERMI(X)
      IMPLICIT NONE
      DOUBLE PRECISION DELTA, X, DMU, FERMI
      PARAMETER(DMU=10.05D0)
      PARAMETER(DELTA=1.04D0)

      FERMI=1.D0+DEXP((X-DMU)/DELTA)
      FERMI=1.D0/FERMI

      RETURN
      END
C__________________________________________________________________________

      FUNCTION DFERMI(X)
      IMPLICIT NONE
      DOUBLE PRECISION DELTA, X, DMU, ERMI, FERMI, DFERMI
      PARAMETER(DMU=10.05D0)
      PARAMETER(DELTA=1.04D0)

      ERMI=DEXP((X-DMU)/DELTA)
      FERMI=1.D0+ERMI
      DFERMI=-ERMI/DELTA
      DFERMI=DFERMI/(FERMI**2)

      RETURN
      END
C__________________________________________________________________________

      FUNCTION D2FERMI(X)
      IMPLICIT NONE
      DOUBLE PRECISION DELTA, X, DMU, ERMI, FERMI, DFERMI, D2FERMI
      PARAMETER(DMU=10.05D0)
      PARAMETER(DELTA=1.04D0)

      ERMI=DEXP((X-DMU)/DELTA)
      FERMI=1.D0+ERMI
      DFERMI=ERMI*(ERMI-1.D0)/DELTA/DELTA
      D2FERMI=DFERMI/(FERMI**3)

      RETURN
      END
C__________________________________________________________________________

      FUNCTION VDW(X)
      IMPLICIT NONE
      DOUBLE PRECISION C6,C8,C10,C12,X,VDW
      PARAMETER(C6=75600.D0)
      PARAMETER(C8=9122400.D0)
      PARAMETER(C10=2.09D8)
      PARAMETER(C12=7.78D10)

      VDW=-C6/X**6-C8/X**8-C10/X**10-C12/X**12

      RETURN
      END
C__________________________________________________________________________

      FUNCTION DVDW(X)
      IMPLICIT NONE
      DOUBLE PRECISION C6,C8,C10,C12,X,DVDW
      PARAMETER(C6=75600.D0)
      PARAMETER(C8=9122400.D0)
      PARAMETER(C10=2.09D8)
      PARAMETER(C12=7.78D10)

      DVDW=6*C6/X**7+8*C8/X**9+10*C10/X**11+12*C12/X**13

      RETURN
      END
C__________________________________________________________________________

      FUNCTION D2VDW(X)
      IMPLICIT NONE
      DOUBLE PRECISION C6,C8,C10,C12,X,D2VDW
      PARAMETER(C6=75600.D0)
      PARAMETER(C8=9122400.D0)
      PARAMETER(C10=2.09D8)
      PARAMETER(C12=7.78D10)

      D2VDW=-42.D0*C6/X**8-72.D0*C8/X**10
      D2VDW=D2VDW-110.D0*C10/X**12-156.D0*C12/X**14

      RETURN
      END
C__________________________________________________________________________

      FUNCTION XMORSE(X)
      IMPLICIT NONE
      DOUBLE PRECISION DM0,TAU,D0,X,XMORSE
      PARAMETER(DM0=0.3D0)
      PARAMETER(TAU=9.75D0)
      PARAMETER(D0=10.3D0)

      XMORSE=DEXP(TAU*(1.D0-X/D0))
      XMORSE=DM0*XMORSE*(XMORSE-2.D0)

      RETURN
      END
C__________________________________________________________________________

      FUNCTION DMORSE(X)
      IMPLICIT NONE
      DOUBLE PRECISION DM0,TAU,D0,MORSE,X,DMORSE
      PARAMETER(DM0=0.3D0)
      PARAMETER(TAU=9.75D0)
      PARAMETER(D0=10.3D0)

      MORSE=DEXP(TAU*(1.D0-X/D0))
      DMORSE=2.D0*TAU*DM0*MORSE*(1.D0-MORSE)
      DMORSE=DMORSE/D0

      RETURN
      END
C__________________________________________________________________________

      FUNCTION D2MORSE(X)
      IMPLICIT NONE
      DOUBLE PRECISION DM0,TAU,D0,X,D2MORSE
      PARAMETER(DM0=0.3D0)
      PARAMETER(TAU=9.75D0)
      PARAMETER(D0=10.3D0)
      DOUBLE PRECISION MORSE,DMORSE

      MORSE=DEXP(TAU*(1.D0-X/D0))
      DMORSE=2.D0*TAU*TAU*DM0*MORSE*(2.D0*MORSE-1.D0)
      D2MORSE=DMORSE/D0/D0

      RETURN
      END
