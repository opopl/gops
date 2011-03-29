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
      SUBROUTINE AMB(XVEC,GRAD,EREAL,GRADT,SECT)
      USE COMMONS
      USE MODAMBER
      USE MODAMBER2
      USE MODHESS
      IMPLICIT NONE
      DOUBLE PRECISION XVEC(3*NATOMS), GRAD(3*NATOMS), EREAL
      LOGICAL GRADT, SECT
      INTEGER J1, J2
      LOGICAL CONNECTIVITY_KNOWN
      SAVE


      IF (ANG.EQ.0) THEN
        COUNT=1
        CALL AMBERMASS
        CALL AMBERS
      END IF

C     IF (COUNT.EQ.1) THEN
C       DO J1=1,NATOMS
C         XVEC(3*(J1-1)+1)=X(J1)
C         XVEC(3*(J1-1)+2)=Y(J1)
C         XVEC(3*(J1-1)+3)=Z(J1)
C       ENDDO
C     ELSE
        DO J1=1,NATOMS
           X(J1)=XVEC(3*(J1-1)+1)
           Y(J1)=XVEC(3*(J1-1)+2)
           Z(J1)=XVEC(3*(J1-1)+3)
        ENDDO
C     END IF
      COUNT=0

      CALL AMBERENERGY

      EREAL=TOTENERGY

      IF (GRADT) THEN
         CALL AMBG
         DO J1=1,NATOMS
            GRAD(3*(J1-1)+1)=DEBYDX(J1)
            GRAD(3*(J1-1)+2)=DEBYDY(J1)
            GRAD(3*(J1-1)+3)=DEBYDZ(J1)
         ENDDO
      ENDIF

      IF (SECT) THEN
         CALL SECONDDERIVS
         DO J1=1,3*NATOMS
            DO J2=1,3*NATOMS
               HESS(J2,J1)=HELL(J2,J1)
            ENDDO
         ENDDO
      ENDIF

      RETURN
      END

      SUBROUTINE ONE4
      USE MODAMBER
      USE MODAMBER2
      USE COMMONS
      IMPLICIT NONE
      INTEGER J1, J2

      DO J1=1,NATOMS
         DO J2=1,NATOMS
            ONE_FOUR(J2,J1)=0
         ENDDO
      ENDDO
      DO A=1,T
        ONE_FOUR(DA1(A),DA4(A))=1
      ENDDO

C     PRINT *,"ONE-FOUR RELATIONSHIPS"
C     DO A=1,ATOMS
C       DO B=A+1,ATOMS
C         IF (ONE_FOUR(A,B).EQ.1)  PRINT *,A,"-",B
C       ENDDO
C     ENDDO

      RETURN
      END
C
C  SETUP STUFF ONLY NEEDS DOING ONCE.
C
      SUBROUTINE AMBERS
      USE COMMONS
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE
      INTEGER            MARVIN
      INTEGER            CHIRATOM(4),CANINE
      COMMON /CHIR/      CANINE
      LOGICAL CONNECTIVITY_KNOWN

      NDIHEDRALS=0
      INQUIRE (FILE='AMBER_CONNECTIVITY',EXIST=CONNECTIVITY_KNOWN)

      CONNECTIVITY_KNOWN=.FALSE.

      IF (CONNECTIVITY_KNOWN) THEN
         PRINT *,'CONNECTIVITY PREVIOUSLY DETERMINED - READING FILE'
         STOP
      ELSE
C 
C CREATE ENTRY IN ANGARRAY
C 
         DO B=1,ATOMS
            DO C=1,ATOMS
               DO D=B,ATOMS
                  IF (B.NE.C .AND. C.NE.D .AND. B.NE.D) THEN
                     IF (.NOT.(BONDS(B,C).NE.1 .OR. BONDS(C,D).NE.1)) THEN
                        ANG=ANG+1
                        AA1(ANG)=B
                        AA2(ANG)=C
                        AA3(ANG)=D
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
         ENDDO

         DO I=1,ATOMS
            DO J=1,ATOMS
               DO K=1,ATOMS
                  DO L=I,ATOMS
                     IF (BONDS(I,J).NE.1 .OR. BONDS(J,I).NE.1) GOTO 113
                     IF (BONDS(J,K).NE.1 .OR. BONDS(K,J).NE.1) GOTO 113
                     IF (BONDS(K,L).NE.1 .OR. BONDS(L,K).NE.1) GOTO 113
                     IF (I.EQ.J .OR. I.EQ.K .OR. I.EQ.L .OR. J.EQ.K .OR. J.EQ.L .OR. K.EQ.L) GOTO 113
     
                     COLIN=0
                     AX=X(J)-X(I)
                     AY=Y(J)-Y(I)
                     AZ=Z(J)-Z(I)
                     BX=X(K)-X(J)
                     BY=Y(K)-Y(J)
                     BZ=Z(K)-Z(J)
                     CX=X(J)-X(L)
                     CY=Y(J)-Y(L)
                     CZ=Z(J)-Z(L)
                     DX=X(L)-X(K)
                     DY=Y(L)-Y(K)
                     DZ=Z(L)-Z(K)
C 
C CHECKING FOR COLINEARITY,CONNECTIVITY AND COINCIDENCE
C 
                     IF (AX*BY.LT.(BX*AY+TINY) .AND. AX*BY.GT.(BX*AY-TINY) .AND. AY*BZ.LT.(BY*AZ+TINY) 
     1                  .AND. AY*BZ.GT.(BY*AZ-TINY) 
     2                  .AND. AX*BZ.LT.(BX*AZ+TINY) .AND. AX*BZ.GT.(AZ*BX-TINY)) COLIN=1
                     IF (BX*DY.LT.(DX*BY+TINY) .AND. BX*DY.GT.(DX*BY-TINY) .AND. BY*DZ.LT.(DY*BZ+TINY) 
     1                  .AND. BY*DZ.GT.(DY*BZ-TINY) 
     2                  .AND. BX*DZ.LT.(BZ*DX+TINY) .AND. BX*DZ.GT.(BZ*DX+TINY)) COLIN=2
                     IF (COLIN.EQ.1) THEN
                       PRINT *,'THREE SITES',I,J,K,'ARE COLINEAR'
                     ELSE IF (COLIN.EQ.2) THEN
                       PRINT *,'THREE SITES',J,K,L,'ARE COLINEAR'
                     ELSE
C 
C CREATE ENTRY IN TORSARRAY
C
                       T=T+1
                       DA1(T)=I
                       DA2(T)=J
                       DA3(T)=K
                       DA4(T)=L
                     ENDIF
113                  CONTINUE
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         CALL ONE4

         PRINT *,' TORSIONS ASSIGNED'

         DO A=1,T
            B=TYPE(DA1(A))
            C=TYPE(DA2(A))
            D=TYPE(DA3(A))        
            E=TYPE(DA4(A))
            MATCH=0
            DO F=1,60
               IF (GENTORSPARAMS(F,1).LT.0.9) GOTO 5
               G=GENTORSPARAMS(F,1)
               H=GENTORSPARAMS(F,2)
               IF (C.EQ.G .AND. D.EQ.H) THEN
                  MATCH=1
                  DID(A)=GENTORSPARAMS(F,3)
                  DVN(A)=GENTORSPARAMS(F,4)
                  DDELTA(A)=GENTORSPARAMS(F,5)
                  DN(A)=GENTORSPARAMS(F,6)
                  DVN2(A)=0.0
                  DVN3(A)=0.0
                  DDELTA2(A)=0.0
                  DDELTA3(A)=0.0
                  DN2(A)=0.0
                  DN3(A)=0.0
               ELSE IF (D.EQ.G .AND. C.EQ.H) THEN
                  MATCH=1
                  DID(A)=GENTORSPARAMS(F,3)
                  DVN(A)=GENTORSPARAMS(F,4)
                  DDELTA(A)=GENTORSPARAMS(F,5)
                  DN(A)=GENTORSPARAMS(F,6)
                  DVN2(A)=0.0
                  DVN3(A)=0.0
                  DDELTA2(A)=0.0
                  DDELTA3(A)=0.0
                  DN2(A)=0.0
                  DN3(A)=0.0
               ENDIF
            ENDDO
5           CONTINUE

            DO F=1,50
               IF (SPECTORSPARAMS(F,1).LT.0.9) GOTO 6
               G=SPECTORSPARAMS(F,1)
               H=SPECTORSPARAMS(F,2)
               I=SPECTORSPARAMS(F,3)
               J=SPECTORSPARAMS(F,4)
               IF (B.EQ.G .AND. C.EQ.H .AND. D.EQ.I .AND. E.EQ.J) THEN
                  MATCH=1
                  DID(A)=SPECTORSPARAMS(F,5)
                  DVN(A)=SPECTORSPARAMS(F,6)
                  DDELTA(A)=SPECTORSPARAMS(F,7)
                  DN(A)=SPECTORSPARAMS(F,8)
                  DVN2(A)=SPECTORSPARAMS(F,9)
                  DDELTA2(A)=SPECTORSPARAMS(F,10)
                  DN2(A)=SPECTORSPARAMS(F,11)
                  DVN3(A)=SPECTORSPARAMS(F,12)
                  DDELTA3(A)=SPECTORSPARAMS(F,13)
                  DN3(A)=SPECTORSPARAMS(F,14)
               ELSE IF (E.EQ.G .AND. D.EQ.H .AND. C.EQ.I .AND. B.EQ.J) THEN
                  MATCH=1
                  DID(A)=SPECTORSPARAMS(F,5)
                  DVN(A)=SPECTORSPARAMS(F,6)
                  DDELTA(A)=SPECTORSPARAMS(F,7)
                  DN(A)=SPECTORSPARAMS(F,8)
                  DVN2(A)=SPECTORSPARAMS(F,9)
                  DDELTA2(A)=SPECTORSPARAMS(F,10)
                  DN2(A)=SPECTORSPARAMS(F,11)
                  DVN3(A)=SPECTORSPARAMS(F,12)
                  DDELTA3(A)=SPECTORSPARAMS(F,13)
                  DN3(A)=SPECTORSPARAMS(F,14)
               ENDIF
            ENDDO
6           CONTINUE
         ENDDO

         IF (MATCH.EQ.0) THEN
            DID(A)=1
            DVN(A)=0.0
            DVN2(A)=0.0
            DVN3(A)=0.0
         ENDIF

         DO I=1,ATOMS
          DO J=I+1,ATOMS
           DO K=J+1,ATOMS
            DO L=K+1,ATOMS
             IF (.NOT.(I.EQ.J .OR. I.EQ.K .OR. I.EQ.L .OR. J.EQ.K .OR. J.EQ.L .OR. K.EQ.L)) THEN
C 
C THIS IF CONSTRUCT PUTS THE ATOMS IN THE IMPROPER IN THE CORRECT ORDER 1-2-3-4
C 
             IF (BONDS(I,K).EQ.1 .AND. BONDS(J,K).EQ.1 .AND. BONDS(K,L).EQ.1) THEN
              IMP=IMP+1
              IF ((TYPE(I) .LE. TYPE(J)) .AND. (TYPE(J) .LE. TYPE(L))) THEN 
               IA1(IMP)=I
               IA2(IMP)=J
               IA3(IMP)=K
               IA4(IMP)=L
              ELSE IF ((TYPE(I) .LE. TYPE(L)) .AND. (TYPE(L) .LT. TYPE(J))) THEN
               IA1(IMP)=I
               IA2(IMP)=L
               IA3(IMP)=K
               IA4(IMP)=J
              ELSE IF ((TYPE(J) .LT. TYPE(I)) .AND. (TYPE(I) .LE. TYPE(L))) THEN
               IA1(IMP)=J
               IA2(IMP)=I
               IA3(IMP)=K
               IA4(IMP)=L
              ELSE IF ((TYPE(J) .LE. TYPE(L)) .AND. (TYPE(L) .LT. TYPE(I))) THEN
               IA1(IMP)=J
               IA2(IMP)=L
               IA3(IMP)=K
               IA4(IMP)=I
              ELSE IF ((TYPE(L) .LT. TYPE(I)) .AND. (TYPE(I) .LE. TYPE(J))) THEN
               IA1(IMP)=L
               IA2(IMP)=I
               IA3(IMP)=K
               IA4(IMP)=J
              ELSE IF ((TYPE(L) .LT. TYPE(J)) .AND. (TYPE(J) .LT. TYPE(I))) THEN
               IA1(IMP)=L
               IA2(IMP)=J
               IA3(IMP)=K
               IA4(IMP)=I
              ENDIF
            
             ELSE IF (BONDS(I,J).EQ.1 .AND. BONDS(I,K).EQ.1 .AND. BONDS(I,L).EQ.1) THEN
              IMP=IMP+1
              IF ((TYPE(J) .LE. TYPE(K)) .AND. (TYPE(K) .LE. TYPE(L))) THEN 
               IA1(IMP)=J
               IA2(IMP)=K
               IA3(IMP)=I
               IA4(IMP)=L
              ELSE IF ((TYPE(J) .LE. TYPE(L)) .AND. (TYPE(L) .LT. TYPE(K))) THEN
               IA1(IMP)=J
               IA2(IMP)=L
               IA3(IMP)=I
               IA4(IMP)=K
              ELSE IF ((TYPE(K) .LT. TYPE(J)) .AND. (TYPE(J) .LE. TYPE(L))) THEN
               IA1(IMP)=K
               IA2(IMP)=J
               IA3(IMP)=I
               IA4(IMP)=L
              ELSE IF ((TYPE(K) .LE. TYPE(L)) .AND. (TYPE(L) .LT. TYPE(J))) THEN
               IA1(IMP)=K
               IA2(IMP)=L
               IA3(IMP)=I
               IA4(IMP)=J
              ELSE IF ((TYPE(L) .LT. TYPE(J)) .AND. (TYPE(J) .LE. TYPE(K))) THEN
               IA1(IMP)=L
               IA2(IMP)=J
               IA3(IMP)=I
               IA4(IMP)=K
              ELSE IF ((TYPE(L) .LT. TYPE(K)) .AND. (TYPE(K) .LT. TYPE(J))) THEN
               IA1(IMP)=L
               IA2(IMP)=K
               IA3(IMP)=I
               IA4(IMP)=J
              ENDIF
   
             ELSE IF (BONDS(I,J).EQ.1 .AND. BONDS(J,K).EQ.1 .AND. BONDS(J,L).EQ.1) THEN
              IMP=IMP+1
              IF ((TYPE(I) .LE. TYPE(K)) .AND. (TYPE(K) .LE. TYPE(L))) THEN 
               IA1(IMP)=I
               IA2(IMP)=K
               IA3(IMP)=J
               IA4(IMP)=L
              ELSE IF ((TYPE(I) .LE. TYPE(L)) .AND. (TYPE(L) .LT. TYPE(K))) THEN
               IA1(IMP)=I
               IA2(IMP)=L
               IA3(IMP)=J
               IA4(IMP)=K
              ELSE IF ((TYPE(K) .LT. TYPE(I)) .AND. (TYPE(I) .LE. TYPE(L))) THEN
               IA1(IMP)=K
               IA2(IMP)=I
               IA3(IMP)=J
               IA4(IMP)=L
              ELSE IF ((TYPE(K) .LE. TYPE(L)) .AND. (TYPE(L) .LT. TYPE(I))) THEN
               IA1(IMP)=K
               IA2(IMP)=L
               IA3(IMP)=J
               IA4(IMP)=I
              ELSE IF ((TYPE(L) .LT. TYPE(I)) .AND. (TYPE(I) .LE. TYPE(K))) THEN
               IA1(IMP)=L
               IA2(IMP)=I
               IA3(IMP)=J
               IA4(IMP)=K
              ELSE IF ((TYPE(L) .LT. TYPE(K)) .AND. (TYPE(K) .LT. TYPE(I))) THEN
               IA1(IMP)=L
               IA2(IMP)=K
               IA3(IMP)=J
               IA4(IMP)=I
              ENDIF
   
             ELSE IF (BONDS(I,L).EQ.1 .AND. BONDS(J,L).EQ.1 .AND. BONDS(K,L).EQ.1) THEN
              IMP=IMP+1
              IF ((TYPE(I) .LE. TYPE(J)) .AND. (TYPE(J) .LE. TYPE(K))) THEN 
               IA1(IMP)=I
               IA2(IMP)=J
               IA3(IMP)=L
               IA4(IMP)=K
              ELSE IF ((TYPE(I) .LE. TYPE(K)) .AND. (TYPE(K) .LT. TYPE(J))) THEN
               IA1(IMP)=I
               IA2(IMP)=K
               IA3(IMP)=L
               IA4(IMP)=J
              ELSE IF ((TYPE(J) .LT. TYPE(I)) .AND. (TYPE(I) .LE. TYPE(K))) THEN
               IA1(IMP)=J
               IA2(IMP)=I
               IA3(IMP)=L
               IA4(IMP)=K
              ELSE IF ((TYPE(J) .LE. TYPE(K)) .AND. (TYPE(K) .LT. TYPE(I))) THEN
               IA1(IMP)=J
               IA2(IMP)=K
               IA3(IMP)=L
               IA4(IMP)=I
              ELSE IF ((TYPE(K) .LT. TYPE(I)) .AND. (TYPE(I) .LE. TYPE(J))) THEN
               IA1(IMP)=K
               IA2(IMP)=I
               IA3(IMP)=L
               IA4(IMP)=J
              ELSE IF ((TYPE(K) .LT. TYPE(J)) .AND. (TYPE(J) .LT. TYPE(I))) THEN
               IA1(IMP)=K
               IA2(IMP)=J
               IA3(IMP)=L
               IA4(IMP)=I
              ENDIF
              ENDIF
             ENDIF
            ENDDO
           ENDDO
          ENDDO
         ENDDO

         PRINT *,' IMPROPERS ASSIGNED'

C
C ASSIGN PARAMETERS TO EACH IMPROPER
C
         DO A=1,IMP
            B=TYPE(IA1(A))
            C=TYPE(IA2(A))        
            D=TYPE(IA3(A))
            E=TYPE(IA4(A))
            DO F=1,15
              G=GENIMPPARAMS(F,1)
              H=GENIMPPARAMS(F,2)
              IF (D.EQ.G .AND. (E.EQ.H .OR. C.EQ.H .OR. B.EQ.H)) THEN
                IVN(A)=GENIMPPARAMS(F,3)
                IDELTA(A)=GENIMPPARAMS(F,4)
                IN1(A)=GENIMPPARAMS(F,5)
              ENDIF
            ENDDO
            DO F=1,4
              G=MIDIMPPARAMS(F,1)
              H=MIDIMPPARAMS(F,2)
              I=MIDIMPPARAMS(F,3)
              IF (C.EQ.G .AND. D.EQ.H .AND. E.EQ.I) THEN
                IVN(A)=MIDIMPPARAMS(F,4)
                IDELTA(A)=MIDIMPPARAMS(F,5)
                IN1(A)=MIDIMPPARAMS(F,6)
              ELSE IF (B.EQ.G .AND. C.EQ.I .AND. D.EQ.H) THEN
                IVN(A)=MIDIMPPARAMS(F,4)
                IDELTA(A)=MIDIMPPARAMS(F,5)
                IN1(A)=MIDIMPPARAMS(F,6)
              ENDIF
            ENDDO
            DO F=1,15
              G=SPECIMPPARAMS(F,1)
              H=SPECIMPPARAMS(F,2)
              I=SPECIMPPARAMS(F,3)
              J=SPECIMPPARAMS(F,4)
      
              IF (B.EQ.G .AND. C.EQ.H .AND. D.EQ.I .AND. E.EQ.J) THEN
                IVN(A)=SPECIMPPARAMS(F,5)
                IDELTA(A)=SPECIMPPARAMS(F,6)
                IN1(A)=SPECIMPPARAMS(F,7)
              ENDIF
            ENDDO
         ENDDO
      END IF

      DO A=1,ATOMS
         DO B=A+1,ATOMS
            RSTAR=VDWR(TYPE(A))+VDWR(TYPE(B))
            EPSILON=SQRT(VDWE(TYPE(A))*VDWE(TYPE(B)))
            VDWA(A,B)=EPSILON*RSTAR**12
            VDWB(A,B)=2*EPSILON*RSTAR**6  
         ENDDO
      ENDDO

      CANINE=1
      DO A=1,ATOMS
        IF (CHIRAL(A).EQ.1) GOTO 111

        IF (TYPECH(A).NE.'CT') GOTO 112
        MARVIN=1
        DO B=1,ATOMS
          IF (BONDS(A,B).EQ.1 .OR. BONDS(B,A).EQ.1) THEN
            CHIRATOM(MARVIN)=B
            MARVIN=MARVIN+1
          END IF
        END DO

        I=CHIRATOM(1)
        J=CHIRATOM(2)
        K=CHIRATOM(3)
        L=CHIRATOM(4)

        IF (TYPECH(I).NE.TYPECH(J) .AND. TYPECH(I).NE.TYPECH(K) .AND. TYPECH(I).NE.TYPECH(L) .AND. 
     1      TYPECH(J).NE.TYPECH(K) .AND. TYPECH(J).NE.TYPECH(L) .AND. TYPECH(K).NE.TYPECH(L)) THEN
          CHIRAL(A)=1
        ELSE
          GOTO 112
        END IF
111     CONTINUE

        CHIRALARRAY(CANINE,1)=A
        CHIRALARRAY(CANINE,2)=I
        CHIRALARRAY(CANINE,3)=J
        CHIRALARRAY(CANINE,4)=K
        CHIRALARRAY(CANINE,5)=L
        CANINE=CANINE+1
   
112   CONTINUE
      END DO

      PRINT *,'SET UP ROUTINE COMPLETED'

      CANINE=CANINE-1
      DO A=1,CANINE
        I=CHIRALARRAY(A,2)
        J=CHIRALARRAY(A,3)
        K=CHIRALARRAY(A,4)
        L=CHIRALARRAY(A,5)
C        PRINT *,'ATOM',CHIRALARRAY(A,1),' IS CHIRAL AND IS BONDED TO ATOMS',I,' ',J,' ',K,' ',L
      END DO

      DO A=1,ATOMS
        IF (LABEL(A).EQ.'D' .AND. TYPECH(A).NE.'N3') THEN
          DO B=A+1,ATOMS
            IF (LABEL(B).EQ.'D' .AND. BONDS(A,B).EQ.1) THEN
              NDIHEDRALS=NDIHEDRALS+1
              DATOM1(NDIHEDRALS)=A
              DATOM2(NDIHEDRALS)=B
C              PRINT *,'DIHEDRAL',NDIHEDRALS,'  ',A,'-',B
C              PRINT *,'         ',TYPECH(A),' - ',TYPECH(B)
            END IF
          END DO
        END IF
      END DO

      CALL ACONNECTDUMP


      RETURN
      END
C
C  DISTANCES
C  VAN DER WAALS AND CHARGE TERMS
C
      SUBROUTINE AMBERD
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE
      INTEGER NDUMMY1, NDUMMY3
      DOUBLE PRECISION DUMMY, DUMMY2, VIXEN1, VIXEN2, VIXEN3, VIXEN4, VIXEN5

      DO A=1,ATOMS
         VIXEN4=PQ(A)
         DO B=A+1,ATOMS
            DUMMY2=(X(A)-X(B))**2+(Y(A)-Y(B))**2+(Z(A)-Z(B))**2
            VIXEN5=DSQRT(DUMMY2)
            R(A,B)=VIXEN5
            R(B,A)=VIXEN5

            NDUMMY1=1-BONDS(A,B)
            NDUMMY3=NDUMMY1*(1-ONE_THREE(A,B))
            IF (NDUMMY3.NE.0) THEN
               VIXEN3=DBLE(ONE_FOUR(A,B))

               DUMMY=1.0D0/DUMMY2**3
               VIXEN1=(VDWA(A,B)*DUMMY-VDWB(A,B))*DUMMY
               IF (FAKEWATER) THEN
                  VIXEN2=VIXEN4*PQ(B)/(DIELEC*VIXEN5**2)
               ELSE
                  VIXEN2=VIXEN4*PQ(B)/(VIXEN5*DIELEC)
               END IF

               VDWENERGY=VDWENERGY+VIXEN1*(1.0D0-VIXEN3/2.0D0)
               QENERGY=    QENERGY+VIXEN2*(1.0D0-VIXEN3/6.0D0)
            ENDIF

C           VIXEN3=FLOAT(ONE_FOUR(A,B))

C           DUMMY=1.0D0/DUMMY2**3
C           VIXEN1=(VDWA(A,B)*DUMMY-VDWB(A,B))*DUMMY
C           VIXEN2=VIXEN4*PQ(B)/R(A,B)

C           VDWENERGY=VDWENERGY+VIXEN1*DUMMY3*(1.0D0-VIXEN3/2.0D0)
C           QENERGY=    QENERGY+VIXEN2*DUMMY3*(1.0D0-VIXEN3/6.0D0)

         ENDDO
      ENDDO
      QENERGY=QENERGY/DIELEC

      RETURN
      END


      SUBROUTINE CHIRALTEST(CTEST)
      USE COMMONS
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE

      INTEGER  CHIRATOM(4)
      INTEGER  MARVIN
      DOUBLE PRECISION  TEMPXH,TEMPYH,TEMPZH,TEMPXI,TEMPYI,TEMPZI,TEMPXJ,TEMPYJ,TEMPZJ
      DOUBLE PRECISION  TEMPXK,TEMPYK,TEMPZK,TEMPXL,TEMPYL,TEMPZL
      DOUBLE PRECISION  NEWXJ,NEWYJ,NEWZJ,NEWXK,NEWYK,NEWZK,GAMMA,DELTA,DOT
      DOUBLE PRECISION  CROSS(3)
      LOGICAL  CTEST
      INTEGER CANINE
      COMMON /CHIR/ CANINE

      CTEST=.FALSE.

      DO A=1,ATOMS
        IF (CHIRAL(A).NE.1) GOTO 17
        MARVIN=1
        DO B=1,ATOMS
          IF (BONDS(A,B).EQ.1 .OR. BONDS(B,A).EQ.1) THEN
            CHIRATOM(MARVIN)=B
            MARVIN=MARVIN+1
          END IF
        END DO    

        I=CHIRATOM(1)
        J=CHIRATOM(2)
        K=CHIRATOM(3)
        L=CHIRATOM(4)

        IF (TYPECH(I).EQ.'HC' .OR. TYPECH(J).EQ.'HC' .OR. TYPECH(K).EQ.'HC' .OR. TYPECH(L).EQ.'HC') THEN
C
C  WE ARE DEALING WITH THE BETA CARBON OF ISOLEUCINE
C
        ELSE IF (TYPECH(I).EQ.'OH'.OR.TYPECH(J).EQ.'OH'.OR.TYPECH(K).EQ.'OH'.OR.TYPECH(L).EQ.'OH') THEN
C
C  WE ARE DEALING WITH THE BETA CARBON OF THREONINE
C  NOTE THAT WE ORDER THE ATOMS BACKWARDS AS THIS CARBON SHOULD BE 'R' CONFIGURATION
C
          PRINT *,' FOUND THREONINE BETA CARBON'
          IF (TYPECH(I).EQ.'OH') THEN
            CHIRATOM(4)=I
          ELSE IF (TYPECH(J).EQ.'OH') THEN
            CHIRATOM(4)=J
          ELSE IF (TYPECH(K).EQ.'OH') THEN
            CHIRATOM(4)=K
          ELSE
            CHIRATOM(4)=L
          END IF
          IF (TYPECH(I).EQ.'H1') THEN
            CHIRATOM(1)=I
          ELSE IF (TYPECH(J).EQ.'H1') THEN
            CHIRATOM(1)=J
          ELSE IF (TYPECH(K).EQ.'H1') THEN
            CHIRATOM(1)=K
          ELSE
            CHIRATOM(1)=L
          END IF
          IF (LABEL(I).EQ.'D') THEN
            CHIRATOM(3)=I
          ELSE IF (LABEL(J).EQ.'D') THEN
            CHIRATOM(3)=J
          ELSE IF (LABEL(K).EQ.'D') THEN
            CHIRATOM(3)=K
          ELSE 
            CHIRATOM(3)=L
          END IF
          IF ((CHIRATOM(1).EQ.I .OR. CHIRATOM(3).EQ.I .OR. CHIRATOM(4).EQ.I) .AND. (CHIRATOM(1).EQ.J
     1       .OR. CHIRATOM(3).EQ.J .OR. CHIRATOM(4).EQ.J) .AND. (CHIRATOM(1).EQ.K .OR. CHIRATOM(3).EQ.K
     2       .OR. CHIRATOM(4).EQ.K)) THEN
            CHIRATOM(2)=L
          ELSE IF ((CHIRATOM(1).EQ.I .OR. CHIRATOM(3).EQ.I .OR. CHIRATOM(4).EQ.I) .AND. (CHIRATOM(1).EQ.J
     1       .OR. CHIRATOM(3).EQ.J .OR. CHIRATOM(4).EQ.J) .AND. (CHIRATOM(1).EQ.L .OR. CHIRATOM(3).EQ.L
     2       .OR. CHIRATOM(4).EQ.L)) THEN          
            CHIRATOM(2)=K
          ELSE IF ((CHIRATOM(1).EQ.I .OR. CHIRATOM(3).EQ.I .OR. CHIRATOM(4).EQ.I) .AND. (CHIRATOM(1).EQ.K
     1       .OR. CHIRATOM(3).EQ.K .OR. CHIRATOM(4).EQ.K) .AND. (CHIRATOM(1).EQ.L .OR. CHIRATOM(3).EQ.L
     2       .OR. CHIRATOM(4).EQ.L)) THEN          
            CHIRATOM(2)=J
          ELSE 
            CHIRATOM(2)=I
          END IF

          PRINT *,'RANKING FOR THREONINE BETA CARBON'
          PRINT *,'1',TYPECH(CHIRATOM(4)),' ',CHIRATOM(4)
          PRINT *,'2',TYPECH(CHIRATOM(3)),' ',CHIRATOM(3)
          PRINT *,'3',TYPECH(CHIRATOM(2)),' ',CHIRATOM(2)
          PRINT *,'4',TYPECH(CHIRATOM(1)),' ',CHIRATOM(1)

        ELSE
C  WE HAVE AN ALPHA CARBON
          IF (TYPECH(J)(1:1).EQ.'H') THEN
            CHIRATOM(1)=J
          ELSE IF (TYPECH(K)(1:1).EQ.'H') THEN
            CHIRATOM(1)=K
          ELSE IF (TYPECH(L)(1:1).EQ.'H') THEN
            CHIRATOM(1)=L
          END IF
    
          IF ((TYPECH(I).EQ.'N ').OR.(TYPECH(I).EQ.'N3')) THEN
            CHIRATOM(2)=I
          ELSE IF ((TYPECH(K).EQ.'N ').OR.(TYPECH(K).EQ.'N3')) THEN
            CHIRATOM(2)=K
          ELSE IF ((TYPECH(L).EQ.'N ').OR.(TYPECH(L).EQ.'N3')) THEN
            CHIRATOM(2)=L
          END IF
  
          IF (TYPECH(I).EQ.'C ') THEN
            CHIRATOM(3)=I
          ELSE IF (TYPECH(J).EQ.'C ') THEN
            CHIRATOM(3)=J
          ELSE IF (TYPECH(L).EQ.'C ') THEN
            CHIRATOM(3)=L
          END IF
  
          IF (TYPECH(I).EQ.'CT') THEN
            CHIRATOM(4)=I
          ELSE IF (TYPECH(J).EQ.'CT') THEN
            CHIRATOM(4)=J
          ELSE IF (TYPECH(K).EQ.'CT') THEN
            CHIRATOM(4)=K
          END IF
        END IF

          H=A
          I=CHIRATOM(1)
          J=CHIRATOM(2)
          K=CHIRATOM(3)
          L=CHIRATOM(4)
C
C MOVE GROUP TO ORIGIN
C
        TEMPXI=X(I)-X(H)
        TEMPYI=Y(I)-Y(H)
        TEMPZI=Z(I)-Z(H)
        TEMPXJ=X(J)-X(H)
        TEMPYJ=Y(J)-Y(H)
        TEMPZJ=Z(J)-Z(H)
        TEMPXK=X(K)-X(H)
        TEMPYK=Y(K)-Y(H)
        TEMPZK=Z(K)-Z(H)
        TEMPXL=X(L)-X(H)
        TEMPYL=Y(L)-Y(H)
        TEMPZL=Z(L)-Z(H)
        TEMPXH=0.0
        TEMPYH=0.0
        TEMPZH=0.0
  
        GAMMA=-(TEMPXJ*TEMPXI+TEMPYJ*TEMPYI+TEMPZJ*TEMPZI)/(TEMPXI**2+TEMPYI**2+TEMPZI**2)
        DELTA=-(TEMPXK*TEMPXI+TEMPYK*TEMPYI+TEMPZK*TEMPZI)/(TEMPXI**2+TEMPYI**2+TEMPZI**2)       

        NEWXJ=TEMPXJ+GAMMA*TEMPXI
        NEWYJ=TEMPYJ+GAMMA*TEMPYI
        NEWZJ=TEMPZJ+GAMMA*TEMPZI
        NEWXK=TEMPXK+DELTA*TEMPXI
        NEWYK=TEMPYK+DELTA*TEMPYI
        NEWZK=TEMPZK+DELTA*TEMPZI
  
        CROSS(1)=NEWYJ*NEWZK-NEWZJ*NEWYK
        CROSS(2)=NEWZJ*NEWXK-NEWXJ*NEWZK
        CROSS(3)=NEWXJ*NEWYK-NEWYJ*NEWXK
  
        DOT=CROSS(1)*TEMPXI+CROSS(2)*TEMPYI+CROSS(3)*TEMPZI
        IF (DOT .LT. 0.0) THEN
C          PRINT *,'CONFIGURATION IS S'
        ELSE
          CTEST=.TRUE. 
        END IF
17    CONTINUE
      END DO
  
      RETURN
      END

      SUBROUTINE AMBERMASS
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE
      INTEGER M1

      DO M1=1,ATOMS
         IF (TYPECH(M1).EQ.'C ') THEN
            MASS(M1)= 12.01
         ELSE IF (TYPECH(M1).EQ.'CA') THEN
            MASS(M1)= 12.01
         ELSE IF (TYPECH(M1).EQ.'CB') THEN
            MASS(M1)= 12.01
         ELSE IF (TYPECH(M1).EQ.'CC') THEN
            MASS(M1)= 12.01
         ELSE IF (TYPECH(M1).EQ.'CK') THEN
            MASS(M1)= 12.01
         ELSE IF (TYPECH(M1).EQ.'CM') THEN
            MASS(M1)= 12.01
         ELSE IF (TYPECH(M1).EQ.'CN') THEN
            MASS(M1)= 12.01
         ELSE IF (TYPECH(M1).EQ.'CQ') THEN
            MASS(M1)= 12.01
         ELSE IF (TYPECH(M1).EQ.'CR') THEN
            MASS(M1)= 12.01
         ELSE IF (TYPECH(M1).EQ.'CT') THEN
            MASS(M1)= 12.01
         ELSE IF (TYPECH(M1).EQ.'CV') THEN
            MASS(M1)= 12.01
         ELSE IF (TYPECH(M1).EQ.'CW') THEN
            MASS(M1)= 12.01
         ELSE IF (TYPECH(M1).EQ.'C*') THEN
            MASS(M1)= 12.01
         ELSE IF (TYPECH(M1).EQ.'C0') THEN
            MASS(M1)= 12.01
         ELSE IF (TYPECH(M1).EQ.'F ') THEN
            MASS(M1)= 19.00
         ELSE IF (TYPECH(M1).EQ.'H ') THEN
            MASS(M1)= 1.008
         ELSE IF (TYPECH(M1).EQ.'HC') THEN
            MASS(M1)= 1.008
         ELSE IF (TYPECH(M1).EQ.'H1') THEN
            MASS(M1)= 1.008
         ELSE IF (TYPECH(M1).EQ.'H2') THEN
            MASS(M1)= 1.008
         ELSE IF (TYPECH(M1).EQ.'H3') THEN
            MASS(M1)= 1.008
         ELSE IF (TYPECH(M1).EQ.'HA') THEN
            MASS(M1)= 1.008
         ELSE IF (TYPECH(M1).EQ.'H4') THEN
            MASS(M1)= 1.008
         ELSE IF (TYPECH(M1).EQ.'H5') THEN
            MASS(M1)= 1.008
         ELSE IF (TYPECH(M1).EQ.'HO') THEN
            MASS(M1)= 1.008
         ELSE IF (TYPECH(M1).EQ.'HS') THEN
            MASS(M1)= 1.008
         ELSE IF (TYPECH(M1).EQ.'HW') THEN
            MASS(M1)= 1.008
         ELSE IF (TYPECH(M1).EQ.'HP') THEN
            MASS(M1)= 1.008
         ELSE IF (TYPECH(M1).EQ.'N ') THEN
            MASS(M1)= 14.01
         ELSE IF (TYPECH(M1).EQ.'NA') THEN
            MASS(M1)= 14.01
         ELSE IF (TYPECH(M1).EQ.'NB') THEN
            MASS(M1)= 14.01
         ELSE IF (TYPECH(M1).EQ.'NC') THEN
            MASS(M1)= 14.01
         ELSE IF (TYPECH(M1).EQ.'N2') THEN
            MASS(M1)= 14.01
         ELSE IF (TYPECH(M1).EQ.'N3') THEN
            MASS(M1)= 14.01
         ELSE IF (TYPECH(M1).EQ.'N*') THEN
            MASS(M1)= 14.01
         ELSE IF (TYPECH(M1).EQ.'O ') THEN
            MASS(M1)= 16.00
         ELSE IF (TYPECH(M1).EQ.'OW') THEN
            MASS(M1)= 16.00
         ELSE IF (TYPECH(M1).EQ.'OH') THEN
            MASS(M1)= 16.00
         ELSE IF (TYPECH(M1).EQ.'OS') THEN
            MASS(M1)= 16.00
         ELSE IF (TYPECH(M1).EQ.'O2') THEN
            MASS(M1)= 16.00
         ELSE IF (TYPECH(M1).EQ.'P ') THEN
            MASS(M1)= 30.97
         ELSE IF (TYPECH(M1).EQ.'S ') THEN
            MASS(M1)= 32.06
         ELSE IF (TYPECH(M1).EQ.'SH') THEN
            MASS(M1)= 32.06
         END IF
      END DO


      RETURN

      END 



      SUBROUTINE AMBG
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE
      INTEGER J1
      DOUBLE PRECISION DPBYDX,DPBYDY,DPBYDZ,DUBYDX,DUBYDY,DUBYDZ,DVBYDX,DVBYDY,DVBYDZ,U,V
      DOUBLE PRECISION VIXEN1, VIXEN2, VIXEN3,VIXEN4,VIXEN5,VIXEN6,VIXEN7,VIXEN8
C 
C INITIALISE ALL DERIVATIVES
C 
      DO J1=1,ATOMS
         DBONDEBYDX(J1)=0.0D0
         DBONDEBYDY(J1)=0.0D0
         DBONDEBYDZ(J1)=0.0D0
         DANGEBYDX(J1)=0.0D0
         DANGEBYDY(J1)=0.0D0
         DANGEBYDZ(J1)=0.0D0
         DTORSEBYDX(J1)=0.0D0
         DTORSEBYDY(J1)=0.0D0
         DTORSEBYDZ(J1)=0.0D0
         DVDWEBYDX(J1)=0.0D0
         DVDWEBYDY(J1)=0.0D0
         DVDWEBYDZ(J1)=0.0D0
         DQEBYDX(J1)=0.0D0
         DQEBYDY(J1)=0.0D0
         DQEBYDZ(J1)=0.0D0
      ENDDO

C      PRINT *,'CO-ORDINATES ARE:'
C      WRITE (*,FMT='(A,3F11.6)') 'O',X(1),Y(1),Z(1)
C      WRITE (*,FMT='(A,3F11.6)') 'H',X(2),Y(2),Z(2)
C      WRITE (*,FMT='(A,3F11.6)') 'H',X(3),Y(3),Z(3)

       DO A=1,ATOMS
         DO B=A+1,ATOMS
           IF (BONDS(A,B).EQ.1) THEN
  
             DBONDEBYDX(A)=DBONDEBYDX(A)+(2.0D0*KR(TYPE(A),TYPE(B))*(X(A)-X(B))*(1.0D0-(RO(TYPE(A),TYPE(B))/R(A,B))))
             DBONDEBYDY(A)=DBONDEBYDY(A)+(2.0D0*KR(TYPE(A),TYPE(B))*(Y(A)-Y(B))*(1.0D0-(RO(TYPE(A),TYPE(B))/R(A,B))))
             DBONDEBYDZ(A)=DBONDEBYDZ(A)+(2.0D0*KR(TYPE(A),TYPE(B))*(Z(A)-Z(B))*(1.0D0-(RO(TYPE(A),TYPE(B))/R(A,B))))
             DBONDEBYDX(B)=DBONDEBYDX(B)+(2.0D0*KR(TYPE(A),TYPE(B))*(X(B)-X(A))*(1.0D0-(RO(TYPE(A),TYPE(B))/R(A,B))))
             DBONDEBYDY(B)=DBONDEBYDY(B)+(2.0D0*KR(TYPE(A),TYPE(B))*(Y(B)-Y(A))*(1.0D0-(RO(TYPE(A),TYPE(B))/R(A,B))))
             DBONDEBYDZ(B)=DBONDEBYDZ(B)+(2.0D0*KR(TYPE(A),TYPE(B))*(Z(B)-Z(A))*(1.0D0-(RO(TYPE(A),TYPE(B))/R(A,B))))

           ENDIF
         ENDDO
       ENDDO

C  PRINT *,"BOND DERIVATIVES CALCULATED"
C 
C FIRST CALCULATE ENERGY GRADIENT FROM CENTRAL ATOM IN ANGLE
C 
      DO D=1,ANG    
       A=AA1(D)
        B=AA2(D)
        C=AA3(D)
 
       VIXEN1=(X(B)**2+Y(B)**2+Z(B)**2)-(X(B)*X(C)+Y(B)*Y(C)+Z(B)*Z(C))
       U=((X(C)-X(B))*X(A)+(Y(C)-Y(B))*Y(A)+(Z(C)-Z(B))*Z(A))+VIXEN1
       V=R(A,B)*R(B,C)
 
       DUBYDX=2.0D0*X(B)-X(C)-X(A)
       DUBYDY=2.0D0*Y(B)-Y(C)-Y(A)
       DUBYDZ=2.0D0*Z(B)-Z(C)-Z(A)
       DVBYDX=((X(B)-X(A))*R(B,C)**2+(X(B)-X(C))*R(A,B)**2)/(R(A,B)*R(B,C))
       DVBYDY=((Y(B)-Y(A))*R(B,C)**2+(Y(B)-Y(C))*R(A,B)**2)/(R(A,B)*R(B,C))
       DVBYDZ=((Z(B)-Z(A))*R(B,C)**2+(Z(B)-Z(C))*R(A,B)**2)/(R(A,B)*R(B,C))
       DPBYDX=((V*DUBYDX)-(U*DVBYDX))/V**2
       DPBYDY=((V*DUBYDY)-(U*DVBYDY))/V**2
       DPBYDZ=((V*DUBYDZ)-(U*DVBYDZ))/V**2

       E=TYPE(A)
       F=TYPE(B)
       G=TYPE(C)

       DANGEBYDX(B)=DANGEBYDX(B)+2.0D0*(KT(E,F,G)*(THETA(D)-TO(E,F,G))*(-DPBYDX/SIN(THETA(D))))
       DANGEBYDY(B)=DANGEBYDY(B)+2.0D0*(KT(E,F,G)*(THETA(D)-TO(E,F,G))*(-DPBYDY/SIN(THETA(D))))
       DANGEBYDZ(B)=DANGEBYDZ(B)+2.0D0*(KT(E,F,G)*(THETA(D)-TO(E,F,G))*(-DPBYDZ/SIN(THETA(D))))
C 
C NOW FROM END ATOMS
C 
       DUBYDX=X(C)-X(B)
       DUBYDY=Y(C)-Y(B)
       DUBYDZ=Z(C)-Z(B)
       DVBYDX=(X(A)-X(B))*R(B,C)/R(A,B)
       DVBYDY=(Y(A)-Y(B))*R(B,C)/R(A,B)
       DVBYDZ=(Z(A)-Z(B))*R(B,C)/R(A,B)
       DPBYDX=((V*DUBYDX)-(U*DVBYDX))/V**2
       DPBYDY=((V*DUBYDY)-(U*DVBYDY))/V**2
       DPBYDZ=((V*DUBYDZ)-(U*DVBYDZ))/V**2

       DANGEBYDX(A)=DANGEBYDX(A)+2.0D0*(KT(E,F,G)*(THETA(D)-TO(E,F,G))*(-DPBYDX/SIN(THETA(D))))
       DANGEBYDY(A)=DANGEBYDY(A)+2.0D0*(KT(E,F,G)*(THETA(D)-TO(E,F,G))*(-DPBYDY/SIN(THETA(D))))
       DANGEBYDZ(A)=DANGEBYDZ(A)+2.0D0*(KT(E,F,G)*(THETA(D)-TO(E,F,G))*(-DPBYDZ/SIN(THETA(D))))


       DUBYDX=X(A)-X(B)
       DUBYDY=Y(A)-Y(B)
       DUBYDZ=Z(A)-Z(B)
       DVBYDX=(X(C)-X(B))*R(A,B)/R(B,C)
       DVBYDY=(Y(C)-Y(B))*R(A,B)/R(B,C)
       DVBYDZ=(Z(C)-Z(B))*R(A,B)/R(B,C)
       DPBYDX=((V*DUBYDX)-(U*DVBYDX))/V**2
       DPBYDY=((V*DUBYDY)-(U*DVBYDY))/V**2
       DPBYDZ=((V*DUBYDZ)-(U*DVBYDZ))/V**2

       DANGEBYDX(C)=DANGEBYDX(C)+2.0D0*(KT(E,F,G)*(THETA(D)-TO(E,F,G))*(-DPBYDX/SIN(THETA(D))))
       DANGEBYDY(C)=DANGEBYDY(C)+2.0D0*(KT(E,F,G)*(THETA(D)-TO(E,F,G))*(-DPBYDY/SIN(THETA(D))))
       DANGEBYDZ(C)=DANGEBYDZ(C)+2.0D0*(KT(E,F,G)*(THETA(D)-TO(E,F,G))*(-DPBYDZ/SIN(THETA(D))))

      ENDDO

C  PRINT *,"ANGLE DERIVATIVES CALCULATED"

       DO A=1,ATOMS
         DO B=A+1,ATOMS
           IF (BONDS(A,B).NE.1) THEN
             VIXEN1=(-12.0D0*VDWA(A,B)/R(A,B)**13)+(6.0D0*VDWB(A,B))/R(A,B)**7
             IF (FAKEWATER) THEN
               VIXEN2=2.0D0*PQ(A)*PQ(B)/(DIELEC*R(A,B)**3)
             ELSE
               VIXEN2=PQ(A)*PQ(B)/(DIELEC*R(A,B)**2)
             END IF
             VIXEN3=(X(A)-X(B))/R(A,B)
             VIXEN4=(Y(A)-Y(B))/R(A,B)
             VIXEN5=(Z(A)-Z(B))/R(A,B)
             VIXEN6=(X(B)-X(A))/R(A,B)
             VIXEN7=(Y(B)-Y(A))/R(A,B)
             VIXEN8=(Z(B)-Z(A))/R(A,B)
 
             IF ((ONE_FOUR(A,B).EQ.1).AND.(ONE_THREE(A,B).NE.1)) THEN
               DVDWEBYDX(A)=DVDWEBYDX(A)+0.5D0*(VIXEN1*VIXEN3)
               DVDWEBYDY(A)=DVDWEBYDY(A)+0.5D0*(VIXEN1*VIXEN4)
               DVDWEBYDZ(A)=DVDWEBYDZ(A)+0.5D0*(VIXEN1*VIXEN5)
               DVDWEBYDX(B)=DVDWEBYDX(B)+0.5D0*(VIXEN1*VIXEN6)
               DVDWEBYDY(B)=DVDWEBYDY(B)+0.5D0*(VIXEN1*VIXEN7)
               DVDWEBYDZ(B)=DVDWEBYDZ(B)+0.5D0*(VIXEN1*VIXEN8)
               DQEBYDX(A)=DQEBYDX(A)-(1.0D0/1.2D0)*VIXEN2*VIXEN3
               DQEBYDY(A)=DQEBYDY(A)-(1.0D0/1.2D0)*VIXEN2*VIXEN4
               DQEBYDZ(A)=DQEBYDZ(A)-(1.0D0/1.2D0)*VIXEN2*VIXEN5
               DQEBYDX(B)=DQEBYDX(B)-(1.0D0/1.2D0)*VIXEN2*VIXEN6
               DQEBYDY(B)=DQEBYDY(B)-(1.0D0/1.2D0)*VIXEN2*VIXEN7
               DQEBYDZ(B)=DQEBYDZ(B)-(1.0D0/1.2D0)*VIXEN2*VIXEN8
             ELSE IF (ONE_THREE(A,B).EQ.1) THEN
C              DO NOTHING
             ELSE 
               DVDWEBYDX(A)=DVDWEBYDX(A)+(VIXEN1*VIXEN3)
               DVDWEBYDY(A)=DVDWEBYDY(A)+(VIXEN1*VIXEN4)
               DVDWEBYDZ(A)=DVDWEBYDZ(A)+(VIXEN1*VIXEN5)
               DVDWEBYDX(B)=DVDWEBYDX(B)+(VIXEN1*VIXEN6)
               DVDWEBYDY(B)=DVDWEBYDY(B)+(VIXEN1*VIXEN7)
               DVDWEBYDZ(B)=DVDWEBYDZ(B)+(VIXEN1*VIXEN8)
               DQEBYDX(A)=DQEBYDX(A)-(VIXEN2*VIXEN3)
               DQEBYDY(A)=DQEBYDY(A)-(VIXEN2*VIXEN4)
               DQEBYDZ(A)=DQEBYDZ(A)-(VIXEN2*VIXEN5)
               DQEBYDX(B)=DQEBYDX(B)-(VIXEN2*VIXEN6)
               DQEBYDY(B)=DQEBYDY(B)-(VIXEN2*VIXEN7)
               DQEBYDZ(B)=DQEBYDZ(B)-(VIXEN2*VIXEN8)
             END IF
           END IF
         ENDDO
       ENDDO


      CALL PDERIVS

      CALL IMPDERIVS
 
      DO A=1,ATOMS
       DEBYDX(A)=DBONDEBYDX(A)+DANGEBYDX(A)+DVDWEBYDX(A)+DTORSEBYDX(A)+DIMPEBYDX(A)+DQEBYDX(A)
       DEBYDY(A)=DBONDEBYDY(A)+DANGEBYDY(A)+DVDWEBYDY(A)+DTORSEBYDY(A)+DIMPEBYDY(A)+DQEBYDY(A)
       DEBYDZ(A)=DBONDEBYDZ(A)+DANGEBYDZ(A)+DVDWEBYDZ(A)+DTORSEBYDZ(A)+DIMPEBYDZ(A)+DQEBYDZ(A)
      ENDDO

      IF (COUNT.EQ.1) PRINT *,"DERIVATIVES CALCULATED"
      COUNT=0

      RETURN
      END

      SUBROUTINE PDERIVS
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE
      DOUBLE PRECISION XE,XF,XG,YE,YF,YG,ZE,ZF,ZG,U,V,UP,VP,P
      DOUBLE PRECISION DLAMBDABYDX,DLAMBDABYDY,DLAMBDABYDZ,DMUBYDX,DMUBYDY,DMUBYDZ
      DOUBLE PRECISION DPBYDX,DPBYDY,DPBYDZ,DUBYDX,DUBYDY,DUBYDZ,DVBYDX,DVBYDY,DVBYDZ
      DOUBLE PRECISION DPHIBYDX,DPHIBYDY,DPHIBYDZ,DUPBYDX,DUPBYDY,DUPBYDZ,DVPBYDX,DVPBYDY,DVPBYDZ
      DOUBLE PRECISION VIXEN1, VIXEN2, VIXEN3
      INTEGER J1
C 
C INITIALISE DERIVATIVES
C 
      DO J1=1,ATOMS
         DTORSEBYDX(J1)=0.0D0
         DTORSEBYDY(J1)=0.0D0
         DTORSEBYDZ(J1)=0.0D0
      ENDDO
C 
C CALCULATE DERIVATIVES FOR ALL ATOMS IN POSITION "A"
C 
      DO J=1,T
C 
C INITIALISE INTERMEDIATES
C 
       DUBYDX=0.0D0
       DUBYDY=0.0D0
       DUBYDZ=0.0D0
       DVBYDX=0.0D0
       DVBYDY=0.0D0
       DVBYDZ=0.0D0
       DPBYDX=0.0D0
       DPBYDY=0.0D0
       DUPBYDX=0.0D0
       DUPBYDY=0.0D0
       DUPBYDZ=0.0D0
       DVPBYDX=0.0D0
       DVPBYDY=0.0D0
       DVPBYDZ=0.0D0
       DPBYDZ=0.0D0

       I=DA1(J)
       IF (I.EQ.0) GOTO 20
       A=DA1(J)
       B=DA2(J)
       C=DA3(J)
       D=DA4(J)
C 
C FIRST DEFINE ALL FUNCTIONS AS ON PAPER
C 
       XE=X(C)-X(B)
       YE=Y(C)-Y(B)
       ZE=Z(C)-Z(B)
       LAMBDA=(((X(B)*XE)+(Y(B)*YE)+(Z(B)*ZE))-((X(A)*XE)+(Y(A)*YE)+(Z(A)*ZE)))/(XE**2+YE**2+ZE**2)
       MU=(((X(B)*XE)+(Y(B)*YE)+(Z(B)*ZE))-((X(D)*XE)+(Y(D)*YE)+(Z(D)*ZE)))/(XE**2+YE**2+ZE**2)
       XF=X(B)-X(A)-(LAMBDA*XE)
       YF=Y(B)-Y(A)-(LAMBDA*YE)
       ZF=Z(B)-Z(A)-(LAMBDA*ZE)
       XG=X(B)-X(D)-(MU*XE)
       YG=Y(B)-Y(D)-(MU*YE)
       ZG=Z(B)-Z(D)-(MU*ZE)
       VIXEN1=(X(B)**2+Y(B)**2+Z(B)**2)+MU*((X(A)-X(B))*XE+(Y(A)-Y(B))*YE+(Z(A)-Z(B))*ZE)
       VIXEN2=LAMBDA*((X(D)-X(B))*XE+(Y(D)-Y(B))*YE+(Z(D)-Z(B))*ZE)-((X(B)*X(D))+(Y(B)*Y(D))+(Z(B)*Z(D)))
       VIXEN3=(X(D)-X(B))*X(A)+(Y(D)-Y(B))*Y(A)+(Z(D)-Z(B))*Z(A)
       U=VIXEN1+VIXEN2+VIXEN3+(LAMBDA*MU*R(B,C)**2)
       UP=SQRT(XF**2+YF**2+ZF**2)
       VP=SQRT(XG**2+YG**2+ZG**2)
       V=UP*VP
       P=U/V

       IF (P.LT.-1) P=P+TINY
       IF (P.GT.1) P=P-TINY
C 
C NOW CALCULATE DERIVATIVES
C 
       DLAMBDABYDX=-XE/R(B,C)**2
       DLAMBDABYDY=-YE/R(B,C)**2
       DLAMBDABYDZ=-ZE/R(B,C)**2
       DUBYDX=(MU*XE)+X(D)-X(B)+DLAMBDABYDX*((X(D)-X(B))*XE+(Y(D)-Y(B))*YE+(Z(D)-Z(B))*ZE+(MU*R(B,C)**2))
       DUBYDY=(MU*YE)+Y(D)-Y(B)+DLAMBDABYDY*((X(D)-X(B))*XE+(Y(D)-Y(B))*YE+(Z(D)-Z(B))*ZE+(MU*R(B,C)**2))
       DUBYDZ=(MU*ZE)+Z(D)-Z(B)+DLAMBDABYDZ*((X(D)-X(B))*XE+(Y(D)-Y(B))*YE+(Z(D)-Z(B))*ZE+(MU*R(B,C)**2))

       DVBYDX=(VP/UP)*(-XF-DLAMBDABYDX*((XF*XE)+(YF*YE)+(ZF*ZE)))
       DVBYDY=(VP/UP)*(-YF-DLAMBDABYDY*((XF*XE)+(YF*YE)+(ZF*ZE)))
       DVBYDZ=(VP/UP)*(-ZF-DLAMBDABYDZ*((XF*XE)+(YF*YE)+(ZF*ZE)))

       DPBYDX=((V*DUBYDX)-(U*DVBYDX))/V**2
       DPBYDY=((V*DUBYDY)-(U*DVBYDY))/V**2
       DPBYDZ=((V*DUBYDZ)-(U*DVBYDZ))/V**2
C 
C NOW ADD TERMS TO ENERGY GRADIENTS
C 
       CALL HAIRY(DPBYDX,DPBYDY,DPBYDZ)
C 
C NOW FOR ALL ATOMS IN POSITION "B"
C 
       I=DA2(J)      

       VIXEN1=(R(B,C)**2)*(X(A)+X(C)-2.0D0*X(B))-2.0D0*((X(B)-X(A))*XE+(Y(B)-Y(A))*YE+(Z(B)-Z(A))*ZE)*(X(B)-X(C))
       DLAMBDABYDX=VIXEN1/(R(B,C)**4)
       VIXEN1=(R(B,C)**2)*(Y(A)+Y(C)-2.0D0*Y(B))-2.0D0*((X(B)-X(A))*XE+(Y(B)-Y(A))*YE+(Z(B)-Z(A))*ZE)*(Y(B)-Y(C))
       DLAMBDABYDY=VIXEN1/(R(B,C)**4)
       VIXEN1=(R(B,C)**2)*(Z(A)+Z(C)-2.0D0*Z(B))-2.0D0*((X(B)-X(A))*XE+(Y(B)-Y(A))*YE+(Z(B)-Z(A))*ZE)*(Z(B)-Z(C))
       DLAMBDABYDZ=VIXEN1/(R(B,C)**4)

       VIXEN1=(R(B,C)**2)*(X(D)+X(C)-2.0D0*X(B))-2.0D0*((X(B)-X(D))*XE+(Y(B)-Y(D))*YE+(Z(B)-Z(D))*ZE)*(X(B)-X(C))
       DMUBYDX=VIXEN1/(R(B,C)**4)
       VIXEN1=(R(B,C)**2)*(Y(D)+Y(C)-2.0D0*Y(B))-2.0D0*((X(B)-X(D))*XE+(Y(B)-Y(D))*YE+(Z(B)-Z(D))*ZE)*(Y(B)-Y(C))
       DMUBYDY=VIXEN1/(R(B,C)**4)
       VIXEN1=(R(B,C)**2)*(Z(D)+Z(C)-2.0D0*Z(B))-2.0D0*((X(B)-X(D))*XE+(Y(B)-Y(D))*YE+(Z(B)-Z(D))*ZE)*(Z(B)-Z(C))
       DMUBYDZ=VIXEN1/(R(B,C)**4)

       VIXEN1=MU*(2.0D0*X(B)-X(C)-X(A))+DMUBYDX*((X(A)-X(B))*XE+(Y(A)-Y(B))*YE+(Z(A)-Z(B))*ZE)
       VIXEN2=LAMBDA*(2.0D0*X(B)-X(C)-X(D))+DLAMBDABYDX*((X(D)-X(B))*XE+(Y(D)-Y(B))*YE+(Z(D)-Z(B))*ZE)
       VIXEN3=((LAMBDA*DMUBYDX)+(MU*DLAMBDABYDX))*R(B,C)**2
       DUBYDX=2.0D0*X(B)-X(A)-X(D)+VIXEN1+VIXEN2+VIXEN3+2.0D0*LAMBDA*MU*(X(B)-X(C))
  
       VIXEN1=MU*(2.0D0*Y(B)-Y(C)-Y(A))+DMUBYDY*((X(A)-X(B))*XE+(Y(A)-Y(B))*YE+(Z(A)-Z(B))*ZE)
       VIXEN2=LAMBDA*(2.0D0*Y(B)-Y(C)-Y(D))+DLAMBDABYDY*((X(D)-X(B))*XE+(Y(D)-Y(B))*YE+(Z(D)-Z(B))*ZE)
       VIXEN3=((LAMBDA*DMUBYDY)+(MU*DLAMBDABYDY))*R(B,C)**2
       DUBYDY=2.0D0*Y(B)-Y(A)-Y(D)+VIXEN1+VIXEN2+VIXEN3+2.0D0*LAMBDA*MU*(Y(B)-Y(C))

       VIXEN1=MU*(2.0D0*Z(B)-Z(C)-Z(A))+DMUBYDZ*((X(A)-X(B))*XE+(Y(A)-Y(B))*YE+(Z(A)-Z(B))*ZE)
       VIXEN2=LAMBDA*(2.0D0*Z(B)-Z(C)-Z(D))+DLAMBDABYDZ*((X(D)-X(B))*XE+(Y(D)-Y(B))*YE+(Z(D)-Z(B))*ZE)
       VIXEN3=((LAMBDA*DMUBYDZ)+(MU*DLAMBDABYDZ))*R(B,C)**2
       DUBYDZ=2.0D0*Z(B)-Z(A)-Z(D)+VIXEN1+VIXEN2+VIXEN3+2.0D0*LAMBDA*MU*(Z(B)-Z(C))
       
       DUPBYDX=(XF*(LAMBDA+1)-DLAMBDABYDX*(XE*XF+YE*YF+ZE*ZF))/UP
       DUPBYDY=(YF*(LAMBDA+1)-DLAMBDABYDY*(XE*XF+YE*YF+ZE*ZF))/UP
       DUPBYDZ=(ZF*(LAMBDA+1)-DLAMBDABYDZ*(XE*XF+YE*YF+ZE*ZF))/UP

       DVPBYDX=(XG*(MU+1)-DMUBYDX*(XE*XG+YE*YG+ZE*ZG))/VP
       DVPBYDY=(YG*(MU+1)-DMUBYDY*(XE*XG+YE*YG+ZE*ZG))/VP
       DVPBYDZ=(ZG*(MU+1)-DMUBYDZ*(XE*XG+YE*YG+ZE*ZG))/VP

       DVBYDX=(UP*DVPBYDX)+(VP*DUPBYDX)
       DVBYDY=(UP*DVPBYDY)+(VP*DUPBYDY)
       DVBYDZ=(UP*DVPBYDZ)+(VP*DUPBYDZ)

       DPBYDX=((V*DUBYDX)-(U*DVBYDX))/V**2
       DPBYDY=((V*DUBYDY)-(U*DVBYDY))/V**2
       DPBYDZ=((V*DUBYDZ)-(U*DVBYDZ))/V**2
      
       CALL HAIRY(DPBYDX,DPBYDY,DPBYDZ)
C 
C NOW FOR ALL ATOMS IN POSITION 3
C 
       I=DA3(J)

       DLAMBDABYDX=(R(B,C)**2*(X(B)-X(A))-2.0D0*((X(B)-X(A))*XE+(Y(B)-Y(A))*YE+(Z(B)-Z(A))*ZE)*
     1             (X(C)-X(B)))/R(B,C)**4    
       DLAMBDABYDY=(R(B,C)**2*(Y(B)-Y(A))-2.0D0*((X(B)-X(A))*XE+(Y(B)-Y(A))*YE+(Z(B)-Z(A))*ZE)*
     1             (Y(C)-Y(B)))/R(B,C)**4
       DLAMBDABYDZ=(R(B,C)**2*(Z(B)-Z(A))-2.0D0*((X(B)-X(A))*XE+(Y(B)-Y(A))*YE+(Z(B)-Z(A))*ZE)*
     1             (Z(C)-Z(B)))/R(B,C)**4     

       DMUBYDX=(R(B,C)**2*(X(B)-X(D))-2.0D0*((X(B)-X(D))*XE+(Y(B)-Y(D))*YE+(Z(B)-Z(D))*ZE)*(X(C)-X(B)))
     1         /R(B,C)**4
       DMUBYDY=(R(B,C)**2*(Y(B)-Y(D))-2.0D0*((X(B)-X(D))*XE+(Y(B)-Y(D))*YE+(Z(B)-Z(D))*ZE)*(Y(C)-Y(B)))
     1         /R(B,C)**4
       DMUBYDZ=(R(B,C)**2*(Z(B)-Z(D))-2.0D0*((X(B)-X(D))*XE+(Y(B)-Y(D))*YE+(Z(B)-Z(D))*ZE)*(Z(C)-Z(B)))
     1         /R(B,C)**4

       VIXEN1=MU*(X(A)-X(B))+DMUBYDX*((X(A)-X(B))*XE+(Y(A)-Y(B))*YE+(Z(A)-Z(B))*ZE)
       VIXEN2=((LAMBDA*DMUBYDX)+(MU*DLAMBDABYDX))*R(B,C)**2+2.0D0*LAMBDA*MU*XE
       DUBYDX=VIXEN1+LAMBDA*(X(D)-X(B))+DLAMBDABYDX*((X(D)-X(B))*XE+(Y(D)-Y(B))*YE+(Z(D)-Z(B))*ZE)+VIXEN2

       VIXEN1=MU*(Y(A)-Y(B))+DMUBYDY*((X(A)-X(B))*XE+(Y(A)-Y(B))*YE+(Z(A)-Z(B))*ZE)
       VIXEN2=((LAMBDA*DMUBYDY)+(MU*DLAMBDABYDY))*R(B,C)**2+2.0D0*LAMBDA*MU*YE
       DUBYDY=VIXEN1+LAMBDA*(Y(D)-Y(B))+DLAMBDABYDY*((X(D)-X(B))*XE+(Y(D)-Y(B))*YE+(Z(D)-Z(B))*ZE)+VIXEN2

       VIXEN1=MU*(Z(A)-Z(B))+DMUBYDZ*((X(A)-X(B))*XE+(Y(A)-Y(B))*YE+(Z(A)-Z(B))*ZE)
       VIXEN2=((LAMBDA*DMUBYDZ)+(MU*DLAMBDABYDZ))*R(B,C)**2+2.0D0*LAMBDA*MU*ZE
       DUBYDZ=VIXEN1+LAMBDA*(Z(D)-Z(B))+DLAMBDABYDZ*((X(D)-X(B))*XE+(Y(D)-Y(B))*YE+(Z(D)-Z(B))*ZE)+VIXEN2

       DUPBYDX=(-LAMBDA*XF-DLAMBDABYDX*(XE*XF+YE*YF+ZE*ZF))/UP
       DUPBYDY=(-LAMBDA*YF-DLAMBDABYDY*(XE*XF+YE*YF+ZE*ZF))/UP
       DUPBYDZ=(-LAMBDA*ZF-DLAMBDABYDZ*(XE*XF+YE*YF+ZE*ZF))/UP

       DVPBYDX=(-MU*XG-DMUBYDX*(XE*XG+YE*YG+ZE*ZG))/VP
       DVPBYDY=(-MU*YG-DMUBYDY*(XE*XG+YE*YG+ZE*ZG))/VP
       DVPBYDZ=(-MU*ZG-DMUBYDZ*(XE*XG+YE*YG+ZE*ZG))/VP

       DVBYDX=(UP*DVPBYDX)+(VP*DUPBYDX)
       DVBYDY=(UP*DVPBYDY)+(VP*DUPBYDY)
       DVBYDZ=(UP*DVPBYDZ)+(VP*DUPBYDZ)

       DPBYDX=((V*DUBYDX)-(U*DVBYDX))/V**2
       DPBYDY=((V*DUBYDY)-(U*DVBYDY))/V**2
       DPBYDZ=((V*DUBYDZ)-(U*DVBYDZ))/V**2

       CALL HAIRY(DPBYDX,DPBYDY,DPBYDZ)       
C 
C NOW FOR ALL ATOMS IN POSITION 4
C 
       I=DA4(J)

       DMUBYDX=-XE/R(B,C)**2
       DMUBYDY=-YE/R(B,C)**2
       DMUBYDZ=-ZE/R(B,C)**2

       VIXEN1=LAMBDA*DMUBYDX*R(B,C)**2
       VIXEN2=LAMBDA*DMUBYDY*R(B,C)**2
       VIXEN3=LAMBDA*DMUBYDZ*R(B,C)**2
       DUBYDX=(LAMBDA*XE)+(((X(A)-X(B))*XE+(Y(A)-Y(B))*YE+(Z(A)-Z(B))*ZE)*DMUBYDX)+X(A)-X(B)+VIXEN1
       DUBYDY=(LAMBDA*YE)+(((X(A)-X(B))*XE+(Y(A)-Y(B))*YE+(Z(A)-Z(B))*ZE)*DMUBYDY)+Y(A)-Y(B)+VIXEN2
       DUBYDZ=(LAMBDA*ZE)+(((X(A)-X(B))*XE+(Y(A)-Y(B))*YE+(Z(A)-Z(B))*ZE)*DMUBYDZ)+Z(A)-Z(B)+VIXEN3

       DVBYDX=UP*(-XG-DMUBYDX*(XE*XG+YE*YG+ZE*ZG))/VP
       DVBYDY=UP*(-YG-DMUBYDY*(XE*XG+YE*YG+ZE*ZG))/VP
       DVBYDZ=UP*(-ZG-DMUBYDZ*(XE*XG+YE*YG+ZE*ZG))/VP

       DPBYDX=((V*DUBYDX)-(U*DVBYDX))/V**2
       DPBYDY=((V*DUBYDY)-(U*DVBYDY))/V**2
       DPBYDZ=((V*DUBYDZ)-(U*DVBYDZ))/V**2

       CALL HAIRY(DPBYDX,DPBYDY,DPBYDZ)
      ENDDO 
20    CONTINUE

      RETURN
      END

      SUBROUTINE HAIRY(DPBYDX,DPBYDY,DPBYDZ)
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE
      DOUBLE PRECISION DPBYDX,DPBYDY,DPBYDZ
      DOUBLE PRECISION VIXEN1,VIXEN2,VIXEN3

       E=TYPE(A)
       F=TYPE(B)
       G=TYPE(C)
       H=TYPE(D)
       PK=DVN(J)
       PK2=DVN2(J)
       PK3=DVN3(J)
       PN=DN(J)
       PN2=DN2(J)
       PN3=DN3(J)
       PHASE=DDELTA(J)
       PHASE2=DDELTA2(J)
       PHASE3=DDELTA3(J)
       IDIVF=DID(J)

       VIXEN1=PK*PN/IDIVF
       VIXEN2=PK2*PN2/IDIVF
       VIXEN3=PK3*PN3/IDIVF       

       IF (PN.EQ.2) THEN
         DTORSEBYDX(I)=DTORSEBYDX(I)+VIXEN1*DPBYDX*2.0D0*COS(DPHI(J))*COS(PHASE)
         DTORSEBYDY(I)=DTORSEBYDY(I)+VIXEN1*DPBYDY*2.0D0*COS(DPHI(J))*COS(PHASE)
         DTORSEBYDZ(I)=DTORSEBYDZ(I)+VIXEN1*DPBYDZ*2.0D0*COS(DPHI(J))*COS(PHASE)
       ELSE IF (PN.EQ.3) THEN
         DTORSEBYDX(I)=DTORSEBYDX(I)+VIXEN1*DPBYDX*(3.0D0-4.0D0*(SIN(DPHI(J)))**2)*COS(PHASE)
         DTORSEBYDY(I)=DTORSEBYDY(I)+VIXEN1*DPBYDY*(3.0D0-4.0D0*(SIN(DPHI(J)))**2)*COS(PHASE)
         DTORSEBYDZ(I)=DTORSEBYDZ(I)+VIXEN1*DPBYDZ*(3.0D0-4.0D0*(SIN(DPHI(J)))**2)*COS(PHASE)
       ELSE IF (PN.EQ.4) THEN
         DTORSEBYDX(I)=DTORSEBYDX(I)+VIXEN1*DPBYDX*(8.0D0*COS(DPHI(J))**3-4.0D0*COS(DPHI(J)))*COS(PHASE)
         DTORSEBYDY(I)=DTORSEBYDY(I)+VIXEN1*DPBYDY*(8.0D0*COS(DPHI(J))**3-4.0D0*COS(DPHI(J)))*COS(PHASE)
         DTORSEBYDZ(I)=DTORSEBYDZ(I)+VIXEN1*DPBYDZ*(8.0D0*COS(DPHI(J))**3-4.0D0*COS(DPHI(J)))*COS(PHASE)
       ENDIF

       IF (PN2.EQ.1) THEN
         DTORSEBYDX(I)=DTORSEBYDX(I)+VIXEN2*DPBYDX*COS(PHASE2)
         DTORSEBYDY(I)=DTORSEBYDY(I)+VIXEN2*DPBYDY*COS(PHASE2)
         DTORSEBYDZ(I)=DTORSEBYDZ(I)+VIXEN2*DPBYDZ*COS(PHASE2)
       ELSE IF (PN2.EQ.2) THEN
         DTORSEBYDX(I)=DTORSEBYDX(I)+VIXEN2*DPBYDX*2.0D0*COS(DPHI(J))*COS(PHASE2)
         DTORSEBYDY(I)=DTORSEBYDY(I)+VIXEN2*DPBYDY*2.0D0*COS(DPHI(J))*COS(PHASE2)
         DTORSEBYDZ(I)=DTORSEBYDZ(I)+VIXEN2*DPBYDZ*2.0D0*COS(DPHI(J))*COS(PHASE2)
       ELSE IF (PN2.EQ.3) THEN
         DTORSEBYDX(I)=DTORSEBYDX(I)+VIXEN2*DPBYDX*(3.0D0-4.0D0*(SIN(DPHI(J)))**2)*COS(PHASE2)
         DTORSEBYDY(I)=DTORSEBYDY(I)+VIXEN2*DPBYDY*(3.0D0-4.0D0*(SIN(DPHI(J)))**2)*COS(PHASE2)
         DTORSEBYDZ(I)=DTORSEBYDZ(I)+VIXEN2*DPBYDZ*(3.0D0-4.0D0*(SIN(DPHI(J)))**2)*COS(PHASE2)
       ELSE IF (PN2.EQ.4) THEN
         DTORSEBYDX(I)=DTORSEBYDX(I)+VIXEN2*DPBYDX*(8.0D0*COS(DPHI(J))**3-4.0D0*COS(DPHI(J)))*COS(PHASE2)
         DTORSEBYDY(I)=DTORSEBYDY(I)+VIXEN2*DPBYDY*(8.0D0*COS(DPHI(J))**3-4.0D0*COS(DPHI(J)))*COS(PHASE2)
         DTORSEBYDZ(I)=DTORSEBYDZ(I)+VIXEN2*DPBYDZ*(8.0D0*COS(DPHI(J))**3-4.0D0*COS(DPHI(J)))*COS(PHASE2)
       ENDIF

       IF (PN3.EQ.1) THEN
         DTORSEBYDX(I)=DTORSEBYDX(I)+VIXEN3*DPBYDX*COS(PHASE3)
         DTORSEBYDY(I)=DTORSEBYDY(I)+VIXEN3*DPBYDY*COS(PHASE3)
         DTORSEBYDZ(I)=DTORSEBYDZ(I)+VIXEN3*DPBYDZ*COS(PHASE3)
       ELSE IF (PN3.EQ.2) THEN
         DTORSEBYDX(I)=DTORSEBYDX(I)+VIXEN3*DPBYDX*2.0D0*COS(DPHI(J))*COS(PHASE3)
         DTORSEBYDY(I)=DTORSEBYDY(I)+VIXEN3*DPBYDY*2.0D0*COS(DPHI(J))*COS(PHASE3)
         DTORSEBYDZ(I)=DTORSEBYDZ(I)+VIXEN3*DPBYDZ*2.0D0*COS(DPHI(J))*COS(PHASE3)
       ELSE IF (PN3.EQ.3) THEN
         DTORSEBYDX(I)=DTORSEBYDX(I)+VIXEN3*DPBYDX*(3.0D0-4.0D0*(SIN(DPHI(J)))**2)*COS(PHASE3)
         DTORSEBYDY(I)=DTORSEBYDY(I)+VIXEN3*DPBYDY*(3.0D0-4.0D0*(SIN(DPHI(J)))**2)*COS(PHASE3)
         DTORSEBYDZ(I)=DTORSEBYDZ(I)+VIXEN3*DPBYDZ*(3.0D0-4.0D0*(SIN(DPHI(J)))**2)*COS(PHASE3)
       ELSE IF (PN3.EQ.4) THEN
         DTORSEBYDX(I)=DTORSEBYDX(I)+VIXEN3*DPBYDX*(8.0D0*COS(DPHI(J))**3-4.0D0*COS(DPHI(J)))*COS(PHASE3)
         DTORSEBYDY(I)=DTORSEBYDY(I)+VIXEN3*DPBYDY*(8.0D0*COS(DPHI(J))**3-4.0D0*COS(DPHI(J)))*COS(PHASE3)
         DTORSEBYDZ(I)=DTORSEBYDZ(I)+VIXEN3*DPBYDZ*(8.0D0*COS(DPHI(J))**3-4.0D0*COS(DPHI(J)))*COS(PHASE3)
       ENDIF

       RETURN
       END

      SUBROUTINE HAIRYIMP(DPBYDX,DPBYDY,DPBYDZ)
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE
      DOUBLE PRECISION DPBYDX,DPBYDY,DPBYDZ,VIXEN1

      E=TYPE(A)
      F=TYPE(B)
      G=TYPE(C)
      H=TYPE(D)
      IPK=IVN(J)
      IPHASE=IDELTA(J)
      IPN=IN1(J)

      VIXEN1=IPK*IPN
   
      IF (IPN.EQ.2) THEN
        DIMPEBYDX(I)=DIMPEBYDX(I)+2.0D0*DPBYDX*VIXEN1*COS(IPHI(J))*COS(IPHASE)       
        DIMPEBYDY(I)=DIMPEBYDY(I)+2.0D0*DPBYDY*VIXEN1*COS(IPHI(J))*COS(IPHASE)       
        DIMPEBYDZ(I)=DIMPEBYDZ(I)+2.0D0*DPBYDZ*VIXEN1*COS(IPHI(J))*COS(IPHASE)       
      ELSE IF (IPN.EQ.0) THEN

      ELSE 
        PRINT *,"IPN FOR THIS TORSION IS NOT CATERED FOR IN DERIVATIVES"
        STOP
      ENDIF

      RETURN
      END
       
      SUBROUTINE IMPDERIVS
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE
      DOUBLE PRECISION XE,XF,XG,YE,YF,YG,ZE,ZF,ZG,U,V,UP,VP,P
      DOUBLE PRECISION DLAMBDABYDX,DLAMBDABYDY,DLAMBDABYDZ,DMUBYDX,DMUBYDY,DMUBYDZ
      DOUBLE PRECISION DPBYDX,DPBYDY,DPBYDZ,DUBYDX,DUBYDY,DUBYDZ,DVBYDX,DVBYDY,DVBYDZ
      DOUBLE PRECISION DPHIBYDX,DPHIBYDY,DPHIBYDZ,DUPBYDX,DUPBYDY,DUPBYDZ,DVPBYDX,DVPBYDY,DVPBYDZ
      DOUBLE PRECISION VIXEN1, VIXEN2, VIXEN3
      INTEGER J1
C 
C INITIALISE DERIVATIVES
C 
      DO J1=1,ATOMS
         DIMPEBYDX(J1)=0.0D0
         DIMPEBYDY(J1)=0.0D0
         DIMPEBYDZ(J1)=0.0D0
      ENDDO
C 
C CALCULATE DERIVATIVES FOR ALL ATOMS IN POSITION "A"
C 
      DO J=1,IMP
C 
C INITIALISE INTERMEDIATES
C 
       DUBYDX=0.0D0
       DUBYDY=0.0D0
       DUBYDZ=0.0D0
       DVBYDX=0.0D0
       DVBYDY=0.0D0
       DVBYDZ=0.0D0
       DPBYDX=0.0D0
       DPBYDY=0.0D0
       DUPBYDX=0.0D0
       DUPBYDY=0.0D0
       DUPBYDZ=0.0D0
       DVPBYDX=0.0D0
       DVPBYDY=0.0D0
       DVPBYDZ=0.0D0
       DPBYDZ=0.0D0

       I=IA1(J)
       IF (I.EQ.0) GOTO 10
       A=IA1(J)
       B=IA2(J)
       C=IA3(J)
       D=IA4(J)
C 
C FIRST DEFINE ALL FUNCTIONS AS ON PAPER
C 
       XE=X(C)-X(B)
       YE=Y(C)-Y(B)
       ZE=Z(C)-Z(B)
       LAMBDA=(((X(B)*XE)+(Y(B)*YE)+(Z(B)*ZE))-((X(A)*XE)+(Y(A)*YE)+(Z(A)*ZE)))/(XE**2+YE**2+ZE**2)
       MU=(((X(B)*XE)+(Y(B)*YE)+(Z(B)*ZE))-((X(D)*XE)+(Y(D)*YE)+(Z(D)*ZE)))/(XE**2+YE**2+ZE**2)
       XF=X(B)-X(A)-(LAMBDA*XE)
       YF=Y(B)-Y(A)-(LAMBDA*YE)
       ZF=Z(B)-Z(A)-(LAMBDA*ZE)
       XG=X(B)-X(D)-(MU*XE)
       YG=Y(B)-Y(D)-(MU*YE)
       ZG=Z(B)-Z(D)-(MU*ZE)
       VIXEN1=(X(B)**2+Y(B)**2+Z(B)**2)+MU*((X(A)-X(B))*XE+(Y(A)-Y(B))*YE+(Z(A)-Z(B))*ZE)
       VIXEN2=LAMBDA*((X(D)-X(B))*XE+(Y(D)-Y(B))*YE+(Z(D)-Z(B))*ZE)-((X(B)*X(D))+(Y(B)*Y(D))+(Z(B)*Z(D)))
       VIXEN3=(X(D)-X(B))*X(A)+(Y(D)-Y(B))*Y(A)+(Z(D)-Z(B))*Z(A)
       U=VIXEN1+VIXEN2+VIXEN3+(LAMBDA*MU*R(B,C)**2)
       UP=SQRT(XF**2+YF**2+ZF**2)
       VP=SQRT(XG**2+YG**2+ZG**2)
       V=UP*VP
       P=U/V

        IF (P.LT.-1) P=P+TINY
        IF (P.GT.1) P=P-TINY
C 
C NOW CALCULATE DERIVATIVES
C 
       DLAMBDABYDX=-XE/R(B,C)**2
       DLAMBDABYDY=-YE/R(B,C)**2
       DLAMBDABYDZ=-ZE/R(B,C)**2
       DUBYDX=(MU*XE)+X(D)-X(B)+DLAMBDABYDX*((X(D)-X(B))*XE+(Y(D)-Y(B))*YE+(Z(D)-Z(B))*ZE+(MU*R(B,C)**2))
       DUBYDY=(MU*YE)+Y(D)-Y(B)+DLAMBDABYDY*((X(D)-X(B))*XE+(Y(D)-Y(B))*YE+(Z(D)-Z(B))*ZE+(MU*R(B,C)**2))
       DUBYDZ=(MU*ZE)+Z(D)-Z(B)+DLAMBDABYDZ*((X(D)-X(B))*XE+(Y(D)-Y(B))*YE+(Z(D)-Z(B))*ZE+(MU*R(B,C)**2))

       DVBYDX=(VP/UP)*(-XF-DLAMBDABYDX*((XF*XE)+(YF*YE)+(ZF*ZE)))
       DVBYDY=(VP/UP)*(-YF-DLAMBDABYDY*((XF*XE)+(YF*YE)+(ZF*ZE)))
       DVBYDZ=(VP/UP)*(-ZF-DLAMBDABYDZ*((XF*XE)+(YF*YE)+(ZF*ZE)))

       DPBYDX=((V*DUBYDX)-(U*DVBYDX))/V**2
       DPBYDY=((V*DUBYDY)-(U*DVBYDY))/V**2
       DPBYDZ=((V*DUBYDZ)-(U*DVBYDZ))/V**2
C 
C NOW ADD TERMS TO ENERGY GRADIENTS
C 
       CALL HAIRYIMP(DPBYDX,DPBYDY,DPBYDZ)
C 
C NOW FOR ALL ATOMS IN POSITION "B"
C 
       I=IA2(J)      

       VIXEN1=(R(B,C)**2)*(X(A)+X(C)-2.0D0*X(B))-2.0D0*((X(B)-X(A))*XE+(Y(B)-Y(A))*YE+(Z(B)-Z(A))*ZE)*(X(B)-X(C))
       DLAMBDABYDX=VIXEN1/(R(B,C)**4)
       VIXEN1=(R(B,C)**2)*(Y(A)+Y(C)-2.0D0*Y(B))-2.0D0*((X(B)-X(A))*XE+(Y(B)-Y(A))*YE+(Z(B)-Z(A))*ZE)*(Y(B)-Y(C))
       DLAMBDABYDY=VIXEN1/(R(B,C)**4)
       VIXEN1=(R(B,C)**2)*(Z(A)+Z(C)-2.0D0*Z(B))-2.0D0*((X(B)-X(A))*XE+(Y(B)-Y(A))*YE+(Z(B)-Z(A))*ZE)*(Z(B)-Z(C))
       DLAMBDABYDZ=VIXEN1/(R(B,C)**4)

       VIXEN1=(R(B,C)**2)*(X(D)+X(C)-2.0D0*X(B))-2.0D0*((X(B)-X(D))*XE+(Y(B)-Y(D))*YE+(Z(B)-Z(D))*ZE)*(X(B)-X(C))
       DMUBYDX=VIXEN1/(R(B,C)**4)
       VIXEN1=(R(B,C)**2)*(Y(D)+Y(C)-2.0D0*Y(B))-2.0D0*((X(B)-X(D))*XE+(Y(B)-Y(D))*YE+(Z(B)-Z(D))*ZE)*(Y(B)-Y(C))
       DMUBYDY=VIXEN1/(R(B,C)**4)
       VIXEN1=(R(B,C)**2)*(Z(D)+Z(C)-2.0D0*Z(B))-2.0D0*((X(B)-X(D))*XE+(Y(B)-Y(D))*YE+(Z(B)-Z(D))*ZE)*(Z(B)-Z(C))
       DMUBYDZ=VIXEN1/(R(B,C)**4)

       VIXEN1=MU*(2.0D0*X(B)-X(C)-X(A))+DMUBYDX*((X(A)-X(B))*XE+(Y(A)-Y(B))*YE+(Z(A)-Z(B))*ZE)
       VIXEN2=LAMBDA*(2.0D0*X(B)-X(C)-X(D))+DLAMBDABYDX*((X(D)-X(B))*XE+(Y(D)-Y(B))*YE+(Z(D)-Z(B))*ZE)
       VIXEN3=((LAMBDA*DMUBYDX)+(MU*DLAMBDABYDX))*R(B,C)**2
       DUBYDX=2.0D0*X(B)-X(A)-X(D)+VIXEN1+VIXEN2+VIXEN3+2.0D0*LAMBDA*MU*(X(B)-X(C))
  
       VIXEN1=MU*(2.0D0*Y(B)-Y(C)-Y(A))+DMUBYDY*((X(A)-X(B))*XE+(Y(A)-Y(B))*YE+(Z(A)-Z(B))*ZE)
       VIXEN2=LAMBDA*(2.0D0*Y(B)-Y(C)-Y(D))+DLAMBDABYDY*((X(D)-X(B))*XE+(Y(D)-Y(B))*YE+(Z(D)-Z(B))*ZE)
       VIXEN3=((LAMBDA*DMUBYDY)+(MU*DLAMBDABYDY))*R(B,C)**2
       DUBYDY=2.0D0*Y(B)-Y(A)-Y(D)+VIXEN1+VIXEN2+VIXEN3+2.0D0*LAMBDA*MU*(Y(B)-Y(C))

       VIXEN1=MU*(2.0D0*Z(B)-Z(C)-Z(A))+DMUBYDZ*((X(A)-X(B))*XE+(Y(A)-Y(B))*YE+(Z(A)-Z(B))*ZE)
       VIXEN2=LAMBDA*(2.0D0*Z(B)-Z(C)-Z(D))+DLAMBDABYDZ*((X(D)-X(B))*XE+(Y(D)-Y(B))*YE+(Z(D)-Z(B))*ZE)
       VIXEN3=((LAMBDA*DMUBYDZ)+(MU*DLAMBDABYDZ))*R(B,C)**2
       DUBYDZ=2.0D0*Z(B)-Z(A)-Z(D)+VIXEN1+VIXEN2+VIXEN3+2.0D0*LAMBDA*MU*(Z(B)-Z(C))
       
       DUPBYDX=(XF*(LAMBDA+1)-DLAMBDABYDX*(XE*XF+YE*YF+ZE*ZF))/UP
       DUPBYDY=(YF*(LAMBDA+1)-DLAMBDABYDY*(XE*XF+YE*YF+ZE*ZF))/UP
       DUPBYDZ=(ZF*(LAMBDA+1)-DLAMBDABYDZ*(XE*XF+YE*YF+ZE*ZF))/UP

       DVPBYDX=(XG*(MU+1)-DMUBYDX*(XE*XG+YE*YG+ZE*ZG))/VP
       DVPBYDY=(YG*(MU+1)-DMUBYDY*(XE*XG+YE*YG+ZE*ZG))/VP
       DVPBYDZ=(ZG*(MU+1)-DMUBYDZ*(XE*XG+YE*YG+ZE*ZG))/VP

       DVBYDX=(UP*DVPBYDX)+(VP*DUPBYDX)
       DVBYDY=(UP*DVPBYDY)+(VP*DUPBYDY)
       DVBYDZ=(UP*DVPBYDZ)+(VP*DUPBYDZ)

       DPBYDX=((V*DUBYDX)-(U*DVBYDX))/V**2
       DPBYDY=((V*DUBYDY)-(U*DVBYDY))/V**2
       DPBYDZ=((V*DUBYDZ)-(U*DVBYDZ))/V**2
      
       CALL HAIRYIMP(DPBYDX,DPBYDY,DPBYDZ)
C 
C NOW FOR ALL ATOMS IN POSITION 3
C 
       I=IA3(J)

       DLAMBDABYDX=(R(B,C)**2*(X(B)-X(A))-2.0D0*((X(B)-X(A))*XE+(Y(B)-Y(A))*YE+(Z(B)-Z(A))*ZE)*
     1             (X(C)-X(B)))/R(B,C)**4    
       DLAMBDABYDY=(R(B,C)**2*(Y(B)-Y(A))-2.0D0*((X(B)-X(A))*XE+(Y(B)-Y(A))*YE+(Z(B)-Z(A))*ZE)*
     1             (Y(C)-Y(B)))/R(B,C)**4
       DLAMBDABYDZ=(R(B,C)**2*(Z(B)-Z(A))-2.0D0*((X(B)-X(A))*XE+(Y(B)-Y(A))*YE+(Z(B)-Z(A))*ZE)*
     1             (Z(C)-Z(B)))/R(B,C)**4     

       DMUBYDX=(R(B,C)**2*(X(B)-X(D))-2.0D0*((X(B)-X(D))*XE+(Y(B)-Y(D))*YE+(Z(B)-Z(D))*ZE)*(X(C)-X(B)))
     1         /R(B,C)**4
       DMUBYDY=(R(B,C)**2*(Y(B)-Y(D))-2.0D0*((X(B)-X(D))*XE+(Y(B)-Y(D))*YE+(Z(B)-Z(D))*ZE)*(Y(C)-Y(B)))
     1         /R(B,C)**4
       DMUBYDZ=(R(B,C)**2*(Z(B)-Z(D))-2.0D0*((X(B)-X(D))*XE+(Y(B)-Y(D))*YE+(Z(B)-Z(D))*ZE)*(Z(C)-Z(B)))
     1         /R(B,C)**4

       VIXEN1=MU*(X(A)-X(B))+DMUBYDX*((X(A)-X(B))*XE+(Y(A)-Y(B))*YE+(Z(A)-Z(B))*ZE)
       VIXEN2=((LAMBDA*DMUBYDX)+(MU*DLAMBDABYDX))*R(B,C)**2+2.0D0*LAMBDA*MU*XE
       DUBYDX=VIXEN1+LAMBDA*(X(D)-X(B))+DLAMBDABYDX*((X(D)-X(B))*XE+(Y(D)-Y(B))*YE+(Z(D)-Z(B))*ZE)+VIXEN2

       VIXEN1=MU*(Y(A)-Y(B))+DMUBYDY*((X(A)-X(B))*XE+(Y(A)-Y(B))*YE+(Z(A)-Z(B))*ZE)
       VIXEN2=((LAMBDA*DMUBYDY)+(MU*DLAMBDABYDY))*R(B,C)**2+2.0D0*LAMBDA*MU*YE
       DUBYDY=VIXEN1+LAMBDA*(Y(D)-Y(B))+DLAMBDABYDY*((X(D)-X(B))*XE+(Y(D)-Y(B))*YE+(Z(D)-Z(B))*ZE)+VIXEN2

       VIXEN1=MU*(Z(A)-Z(B))+DMUBYDZ*((X(A)-X(B))*XE+(Y(A)-Y(B))*YE+(Z(A)-Z(B))*ZE)
       VIXEN2=((LAMBDA*DMUBYDZ)+(MU*DLAMBDABYDZ))*R(B,C)**2+2.0D0*LAMBDA*MU*ZE
       DUBYDZ=VIXEN1+LAMBDA*(Z(D)-Z(B))+DLAMBDABYDZ*((X(D)-X(B))*XE+(Y(D)-Y(B))*YE+(Z(D)-Z(B))*ZE)+VIXEN2

       DUPBYDX=(-LAMBDA*XF-DLAMBDABYDX*(XE*XF+YE*YF+ZE*ZF))/UP
       DUPBYDY=(-LAMBDA*YF-DLAMBDABYDY*(XE*XF+YE*YF+ZE*ZF))/UP
       DUPBYDZ=(-LAMBDA*ZF-DLAMBDABYDZ*(XE*XF+YE*YF+ZE*ZF))/UP

       DVPBYDX=(-MU*XG-DMUBYDX*(XE*XG+YE*YG+ZE*ZG))/VP
       DVPBYDY=(-MU*YG-DMUBYDY*(XE*XG+YE*YG+ZE*ZG))/VP
       DVPBYDZ=(-MU*ZG-DMUBYDZ*(XE*XG+YE*YG+ZE*ZG))/VP

       DVBYDX=(UP*DVPBYDX)+(VP*DUPBYDX)
       DVBYDY=(UP*DVPBYDY)+(VP*DUPBYDY)
       DVBYDZ=(UP*DVPBYDZ)+(VP*DUPBYDZ)

       DPBYDX=((V*DUBYDX)-(U*DVBYDX))/V**2
       DPBYDY=((V*DUBYDY)-(U*DVBYDY))/V**2
       DPBYDZ=((V*DUBYDZ)-(U*DVBYDZ))/V**2

       CALL HAIRYIMP(DPBYDX,DPBYDY,DPBYDZ)       
C 
C NOW FOR ALL ATOMS IN POSITION 4
C 
       I=IA4(J)

       DMUBYDX=-XE/R(B,C)**2
       DMUBYDY=-YE/R(B,C)**2
       DMUBYDZ=-ZE/R(B,C)**2

       VIXEN1=LAMBDA*DMUBYDX*R(B,C)**2
       VIXEN2=LAMBDA*DMUBYDY*R(B,C)**2
       VIXEN3=LAMBDA*DMUBYDZ*R(B,C)**2
       DUBYDX=(LAMBDA*XE)+(((X(A)-X(B))*XE+(Y(A)-Y(B))*YE+(Z(A)-Z(B))*ZE)*DMUBYDX)+X(A)-X(B)+VIXEN1
       DUBYDY=(LAMBDA*YE)+(((X(A)-X(B))*XE+(Y(A)-Y(B))*YE+(Z(A)-Z(B))*ZE)*DMUBYDY)+Y(A)-Y(B)+VIXEN2
       DUBYDZ=(LAMBDA*ZE)+(((X(A)-X(B))*XE+(Y(A)-Y(B))*YE+(Z(A)-Z(B))*ZE)*DMUBYDZ)+Z(A)-Z(B)+VIXEN3

       DVBYDX=UP*(-XG-DMUBYDX*(XE*XG+YE*YG+ZE*ZG))/VP
       DVBYDY=UP*(-YG-DMUBYDY*(XE*XG+YE*YG+ZE*ZG))/VP
       DVBYDZ=UP*(-ZG-DMUBYDZ*(XE*XG+YE*YG+ZE*ZG))/VP

       DPBYDX=((V*DUBYDX)-(U*DVBYDX))/V**2
       DPBYDY=((V*DUBYDY)-(U*DVBYDY))/V**2
       DPBYDZ=((V*DUBYDZ)-(U*DVBYDZ))/V**2

       CALL HAIRYIMP(DPBYDX,DPBYDY,DPBYDZ)
      ENDDO 
10    CONTINUE

      RETURN
      END


      SUBROUTINE NUMERGRAD
      USE COMMONS
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE
      INTEGER I1
      DOUBLE PRECISION         NUMERDERIV(3*NATOMS)
      DOUBLE PRECISION         ET

      DO I1=1,3*NATOMS
         NUMERDERIV(I1)= 0.0
      END DO

      DO I1=1,ATOMS
         X(I1)= X(I1)+ 0.00001
         CALL AMBERENERGY
         ET= TOTENERGY
         X(I1)= X(I1)- 0.00002
         CALL AMBERENERGY
         NUMERDERIV(3*I1-2)= (ET-TOTENERGY)/0.00002
         X(I1)= X(I1)+ 0.00001

         Y(I1)= Y(I1)+ 0.00001
         CALL AMBERENERGY
         ET= TOTENERGY
         Y(I1)= Y(I1)- 0.00002
         CALL AMBERENERGY
         NUMERDERIV(3*I1-1)= (ET-TOTENERGY)/0.00002
         Y(I1)= Y(I1)+ 0.00001

         Z(I1)= Z(I1)+ 0.00001
         CALL AMBERENERGY
         ET= TOTENERGY
         Z(I1)= Z(I1)- 0.00002
         CALL AMBERENERGY
         NUMERDERIV(3*I1)= (ET-TOTENERGY)/0.00002
         Z(I1)= Z(I1)+ 0.00001
      END DO

      WRITE (*,FMT='(5X,A,5X,A,5X,A,5X,A,5X,A)') 'NUMERICAL','   *   ','ANALYTIC','   *   ','DIFFERENCE'
      WRITE (*,FMT='(5X,A,5X,A,5X,A,5X,A,5X,A)') '---------','   *   ','--------','   *   ','----------'
      DO I1=1,ATOMS
         WRITE (*,FMT='(F20.10,A,F20.10,A,F20.10)') NUMERDERIV(3*I1-2),'  * ',DEBYDX(I1),'   * ',DEBYDX(I1)-NUMERDERIV(3*I1-2)
         WRITE (*,FMT='(F20.10,A,F20.10,A,F20.10)') NUMERDERIV(3*I1-1),'  * ',DEBYDY(I1),'   * ',DEBYDY(I1)-NUMERDERIV(3*I1-1)
         WRITE (*,FMT='(F20.10,A,F20.10,A,F20.10)') NUMERDERIV(3*I1),'  * ',DEBYDZ(I1),'   * ',DEBYDZ(I1)-NUMERDERIV(3*I1)
      END DO


      RETURN

      END

      SUBROUTINE AMBERA
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE

      DO A=1,ANG
         B=AA1(A)
         C=AA2(A)
         D=AA3(A)

         TOP=(X(B)-X(C))*(X(D)-X(C))+(Y(B)-Y(C))*(Y(D)-Y(C))+(Z(B)-Z(C))*(Z(D)-Z(C))
         BOTTOM=R(B,C)*DSQRT((X(D)-X(C))**2+(Y(D)-Y(C))**2+(Z(D)-Z(C))**2)
         THETA(A)=ACOS(TOP/BOTTOM)
      END DO

      DO A=1,T
         I=DA1(A)
         J=DA2(A)
         K=DA3(A)
         L=DA4(A)
   
         COLIN=0
         AX=X(J)-X(I)
         AY=Y(J)-Y(I)
         AZ=Z(J)-Z(I)
         BX=X(K)-X(J)
         BY=Y(K)-Y(J)
         BZ=Z(K)-Z(J)
         CX=X(J)-X(L)
         CY=Y(J)-Y(L)
         CZ=Z(J)-Z(L)
         DX=X(L)-X(K)
         DY=Y(L)-Y(K)
         DZ=Z(L)-Z(K)
   
         IF (AX*BY.LT.(BX*AY+TINY) .AND. AX*BY.GT.(BX*AY-TINY) .AND. AY*BZ.LT.(BY*AZ+TINY) .AND. AY*BZ.GT.(BY*AZ-TINY) 
     1        .AND. AX*BZ.LT.(BX*AZ+TINY) .AND. AX*BZ.GT.(AZ*BX-TINY)) COLIN=1
         IF (BX*DY.LT.(DX*BY+TINY) .AND. BX*DY.GT.(DX*BY-TINY) .AND. BY*DZ.LT.(DY*BZ+TINY) .AND. BY*DZ.GT.(DY*BZ-TINY) 
     1        .AND. BX*DZ.LT.(BZ*DX+TINY) .AND. BX*DZ.GT.(BZ*DX+TINY)) COLIN=2

         IF (COLIN.EQ.1) THEN
            PRINT *,'THREE SITES',I,J,K,'ARE COLINEAR'
         ELSE IF (COLIN.EQ.2) THEN
            PRINT *,'THREE SITES',J,K,L,'ARE COLINEAR'
         ELSE
            LAMBDA=(AX*BX+AY*BY+AZ*BZ)/(BX**2+BY**2+BZ**2)
            MU=(CX*BX+CY*BY+CZ*BZ)/(BX**2+BY**2+BZ**2)
            QX(1)=X(I)+(LAMBDA*BX)
            QY(1)=Y(I)+(LAMBDA*BY)
            QZ(1)=Z(I)+(LAMBDA*BZ)
            QX(2)=X(J)
            QY(2)=Y(J)
            QZ(2)=Z(J)
            QX(3)=X(L)+(MU*BX)
            QY(3)=Y(L)+(MU*BY)
            QZ(3)=Z(L)+(MU*BZ)
            NUMER=(QX(1)-QX(2))*(QX(3)-QX(2))+(QY(1)-QY(2))*(QY(3)-QY(2))+(QZ(1)-QZ(2))*(QZ(3)-QZ(2))
            DENOM=SQRT(((QX(1)-QX(2))**2+(QY(1)-QY(2))**2+(QZ(1)-QZ(2))**2)
     1          *((QX(3)-QX(2))**2+(QY(3)-QY(2))**2+(QZ(3)-QZ(2))**2))
   
            IF ((NUMER/DENOM).LT. -1) NUMER=NUMER+1D-12
            IF ((NUMER/DENOM).GT.1) NUMER=NUMER-1D-12
  
            DPHI(A)=ACOS(NUMER/DENOM)
         ENDIF
      END DO

      DO A=1,IMP 
         E=IA1(A)
         F=IA2(A)
         G=IA3(A)
         H=IA4(A)    
         COLIN=0
         AX=X(F)-X(E)
         AY=Y(F)-Y(E)
         AZ=Z(F)-Z(E)
         BX=X(G)-X(F)
         BY=Y(G)-Y(F)
         BZ=Z(G)-Z(F)
         CX=X(F)-X(H)
         CY=Y(F)-Y(H)
         CZ=Z(F)-Z(H)
         DX=X(H)-X(G)
         DY=Y(H)-Y(G)
         DZ=Z(H)-Z(G)

         IF (AX*BY.LT.(BX*AY+TINY) .AND. AX*BY.GT.(BX*AY-TINY) .AND. AY*BZ.LT.(BY*AZ+TINY) .AND. AY*BZ.GT.(BY*AZ-TINY) 
     1      .AND. AX*BZ.LT.(BX*AZ+TINY) .AND. AX*BZ.GT.(AZ*BX-TINY)) COLIN=1
         IF (BX*DY.LT.(DX*BY+TINY) .AND. BX*DY.GT.(DX*BY-TINY) .AND. BY*DZ.LT.(DY*BZ+TINY) .AND. BY*DZ.GT.(DY*BZ-TINY) 
     1      .AND. BX*DZ.LT.(BZ*DX+TINY) .AND. BX*DZ.GT.(BZ*DX+TINY)) COLIN=2
         IF (COLIN.EQ.1) THEN
            PRINT *,'THREE SITES',E,F,G,' ARE COLINEAR'
         ELSE IF (COLIN.EQ.2) THEN
            PRINT *,'THREE SITES',F,G,H,' ARE COLINEAR'
         ELSE

            LAMBDA=(AX*BX+AY*BY+AZ*BZ)/(BX**2+BY**2+BZ**2)
            MU=(CX*BX+CY*BY+CZ*BZ)/(BX**2+BY**2+BZ**2)
            QX(1)=X(E)+(LAMBDA*BX)
            QY(1)=Y(E)+(LAMBDA*BY)
            QZ(1)=Z(E)+(LAMBDA*BZ)
            QX(2)=X(F)
            QY(2)=Y(F)
            QZ(2)=Z(F)
            QX(3)=X(H)+(MU*BX)
            QY(3)=Y(H)+(MU*BY)
            QZ(3)=Z(H)+(MU*BZ)
              NUMER=(QX(1)-QX(2))*(QX(3)-QX(2))+(QY(1)-QY(2))*(QY(3)-QY(2))+(QZ(1)-QZ(2))*(QZ(3)-QZ(2))
            DENOM=SQRT(((QX(1)-QX(2))**2+(QY(1)-QY(2))**2+(QZ(1)-QZ(2))**2)
     1          *((QX(3)-QX(2))**2+(QY(3)-QY(2))**2+(QZ(3)-QZ(2))**2))

            IF ((NUMER/DENOM).LT. -1) NUMER=NUMER+1D-12
            IF ((NUMER/DENOM).GT.1) NUMER=NUMER-1D-12

            IPHI(A)=ACOS(NUMER/DENOM)
         ENDIF
      END DO

      RETURN

      END
      SUBROUTINE AMBERENERGY
      USE KEY
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE
      DOUBLE PRECISION VIXEN1, VIXEN2, VIXEN3, VIXEN4, VIXEN5, VIXEN6, E1, E2, E3
      CHARACTER FNAMEF*80

C 
C FIRST INITIALISING COUNT VARIABLES
C 
      BENERGY=0.0
      TENERGY=0.0
      PENERGY=0.0
      VDWENERGY=0.0
      TOTENERGY=0.0
      IMPENERGY=0.0
      QENERGY=0.0
C 
C CALCULATING R,BENERGY,VDWENERGY,QENERGY - NOTE THAT WE CALCULATE ALL R AS WE NEED THEM FOR VDW ENERGY
C 
      CALL AMBERD

      DO I=1,BONDNUMBER
         A=BONDARRAY(I,1)
         B=BONDARRAY(I,2)
         BENERGY=BENERGY+KR(TYPE(A),TYPE(B))*((R(A,B)-RO(TYPE(A),TYPE(B)))**2)
      ENDDO

      CALL AMBERA

      DO A=1,ANG
         B=AA1(A)
         C=AA2(A)
         D=AA3(A)

         VIXEN1=KT(TYPE(B),TYPE(C),TYPE(D))*((THETA(A)-TO(TYPE(B),TYPE(C),TYPE(D)))**2)
         TENERGY=TENERGY+VIXEN1
C        PRINT *,'ANGLE ',B,C,D,TYPECH(B),'-',TYPECH(C),'-',TYPECH(D),' THETA=',THETA(A)*57.29577951,' ENERGY=',VIXEN1
      ENDDO

      DO A=1,T
         I=DA1(A)
         J=DA2(A)
         K=DA3(A)
         L=DA4(A)
         VIXEN1=DVN(A)/DID(A)
         VIXEN2=DN(A)*DPHI(A)-DDELTA(A)
 
         VIXEN3=DVN2(A)/DID(A)
         VIXEN4=DN2(A)*DPHI(A)-DDELTA2(A)

         VIXEN5=DVN3(A)/DID(A)
         VIXEN6=DN3(A)*DPHI(A)-DDELTA3(A)

         PENERGY=PENERGY+(VIXEN1*(1+COS(VIXEN2)))+(VIXEN3*(1+COS(VIXEN4)))+(VIXEN5*(1+COS(VIXEN6)))
         E1=VIXEN1*(1+COS(VIXEN2))
         E2=(VIXEN3*(1+COS(VIXEN4)))
         E3=(VIXEN5*(1+COS(VIXEN6)))

C         PRINT *,'TORSION',I,J,K,L,TYPECH(I),'-',TYPECH(J),'-',TYPECH(K),'-',TYPECH(L),' PHI=',DPHI(A)*57.29577951,
C     1           ' ENERGY=',E1+E2+E3
C         PRINT *,'PK=',DVN(A),' IDIVF=',DID(A),' PN=',DN(A),' PHASE=',DDELTA(A)

      ENDDO

      DO A=1,IMP
         IMPENERGY=IMPENERGY+(IVN(A)*(1+COS(IN1(A)*IPHI(A)-IDELTA(A))))
     
C        PRINT *,'IMPROPER',E,F,G,H,' ENERGY= ',IVN(A)*(1+COS(IN1(A)*IPHI(A)-IDELTA(A)))

      ENDDO

      IF (COUNT.EQ.1) THEN 
         PRINT *,'ANG= ',ANG
         PRINT *,'T= ',T
         PRINT *,'IMP= ',IMP
      ENDIF
C     DO A=1,IMP
C        PRINT *,IA1(A),IA2(A),IA3(A),IA4(A)
C        PRINT *,'PHI= ',IMPROPERS(A)%PHI*57.29577951
C     ENDDO

C     PRINT *,'BOND LENGTHS'
C       DO A=1,ATOMS
C        DO B=A,ATOMS
C        IF (BONDS(A,B).EQ.1) PRINT *,A,'-',B,'    ',R(A,B),'  ',TYPECH(A),' ',TYPECH(B)
C        ENDDO
C       ENDDO

C       PRINT *,'BOND ANGLES'
C         DO A=1,ANG
C          PRINT *,ANGLES(A)%A1,'-',ANGLES(A)%A2,'-',ANGLES(A)%A3,' = ',(ANGLES(A)%THETA)*57.29577951
C         ENDDO

C       PRINT *,'TORSION ANGLES'
C       DO A=1,T
C        PRINT *,DA1(A),TORSIONS(A)%A2,TORSIONS(A)%A3,TORSIONS(A)%A4,'    ',&
C                57.29577951*TORSIONS(A)%PHI
C       ENDDO
    
       PRINT *,'TOTAL BOND STRAIN ENERGY= ',BENERGY
       PRINT *,'TOTAL ANGLE STRAIN ENERGY= ',TENERGY
       PRINT *,'TORSION ANGLE ENERGY= ',PENERGY
       PRINT *,'IMPROPER TORSION ENERGY= ',IMPENERGY
       PRINT *,'VDW ENERGY= ',VDWENERGY
       PRINT *,'QENERGY= ',QENERGY
       TOTENERGY=BENERGY+TENERGY+PENERGY+VDWENERGY+IMPENERGY+QENERGY
       PRINT *,'TOTAL ENERGY= ',TOTENERGY
C      IF (COUNT.EQ.1) PRINT *,'ENERGY CALCULATED'

       IF (TOTENERGY.LT.-1.0D4) THEN  
          WRITE(*,'(A)') 'COLD FUSION DIAGNOSED - QUIT'
          IF (FILTH.EQ.0) THEN
             FNAMEF='POINTS.FINAL'
          ELSE
             WRITE(FNAMEF,'(A)') 'POINTS.FINAL.'//TRIM(ADJUSTL(FILTHSTR))
          ENDIF
C         CALL AMBERDUMP(Q,FNAMEF)
          STOP
       ENDIF

      RETURN
      END

      SUBROUTINE AMOVIEDUMP(FRAME)
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE
      INTEGER FRAME, MOVIECOUNT

      WRITE (27,FMT=*) ATOMS
      WRITE (27,FMT='(A,I4)') 'FRAME NO. ',FRAME
      DO MOVIECOUNT=1,ATOMS
         WRITE (27,FMT='(A,3F20.10)') TYPECH(MOVIECOUNT)(1:1), X(MOVIECOUNT), Y(MOVIECOUNT), Z(MOVIECOUNT)
      END DO

      FRAME=FRAME+1

      RETURN

      END
      SUBROUTINE AREAD
      USE COMMONS
      USE KEY
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE
      INTEGER CANINE
      COMMON /CHIR/ CANINE
      DOUBLE PRECISION VIXEN1
      DOUBLE PRECISION CRAP
      PARAMETER (CRAP=18.2223D0)
      CHARACTER(LEN=17) FILTHFILE

      ALLOCATE(DATOM1(NATOMS),DATOM2(NATOMS))
      ALLOCATE(BONDARRAY(NATOMS*4,2))
      ALLOCATE(X(NATOMS),Y(NATOMS),Z(NATOMS),PQ(NATOMS),MASS(NATOMS))
      ALLOCATE(R(NATOMS,NATOMS),VDWA(NATOMS,NATOMS),VDWB(NATOMS,NATOMS))
      ALLOCATE(AA1(NATOMS*5),AA2(NATOMS*5),AA3(NATOMS*5))
      ALLOCATE(THETA(NATOMS*5))
      ALLOCATE(DA1(NATOMS*5),DA2(NATOMS*5),DA3(NATOMS*5),DA4(NATOMS*5))
      ALLOCATE(DPHI(NATOMS*5),DID(NATOMS*5),DVN(NATOMS*5),DVN2(NATOMS*5))
      ALLOCATE(DVN3(NATOMS*5),DDELTA(NATOMS*5),DDELTA2(NATOMS*5),DDELTA3(NATOMS*5))
      ALLOCATE(DN(NATOMS*5),DN2(NATOMS*5),DN3(NATOMS*5))
      ALLOCATE(IA1(NATOMS*5),IA2(NATOMS*5),IA3(NATOMS*5),IA4(NATOMS*5))
      ALLOCATE(IPHI(NATOMS*5),IVN(NATOMS*5),IDELTA(NATOMS*5),IN1(NATOMS*5))
      ALLOCATE(ATNUM(NATOMS),BONDEDTO(NATOMS),TYPE(NATOMS))
      ALLOCATE(BONDS(NATOMS,NATOMS))
      ALLOCATE(ONE_FOUR(NATOMS,NATOMS),ONE_THREE(NATOMS,NATOMS))
      ALLOCATE(LABEL(NATOMS))
      ALLOCATE(TYPECH(NATOMS))
      ALLOCATE(DBONDEBYDX(NATOMS),DBONDEBYDY(NATOMS),DBONDEBYDZ(NATOMS),DANGEBYDX(NATOMS), 
     &                 DANGEBYDY(NATOMS),DANGEBYDZ(NATOMS), 
     &                 DTORSEBYDX(NATOMS),DTORSEBYDY(NATOMS),DTORSEBYDZ(NATOMS),DVDWEBYDX(NATOMS), 
     &                 DVDWEBYDY(NATOMS),DVDWEBYDZ(NATOMS), 
     &                 DEBYDX(NATOMS),DEBYDY(NATOMS),DEBYDZ(NATOMS),DIMPEBYDX(NATOMS), 
     &                 DIMPEBYDY(NATOMS),DIMPEBYDZ(NATOMS), 
     &                 DQEBYDX(NATOMS),DQEBYDY(NATOMS),DQEBYDZ(NATOMS))
      ALLOCATE(BONDHELL(3*NATOMS,3*NATOMS), ANGLEHELL(3*NATOMS,3*NATOMS), TORSHELL(3*NATOMS,3*NATOMS),
     1         IMPHELL(3*NATOMS,3*NATOMS), QHELL(3*NATOMS,3*NATOMS), VDWHELL(3*NATOMS,3*NATOMS),
     2         HELL(3*NATOMS,3*NATOMS))
      ALLOCATE(CHIRAL(NATOMS),CHIRALARRAY(NATOMS,5))

      IF (FILTH.EQ.0) THEN
         OPEN (UNIT=9,FILE='COORDS.AMBER',STATUS='OLD')
      ELSE 
         WRITE (FILTHFILE,'(A)') 'COORDS.AMBER.'//TRIM(ADJUSTL(FILTHSTR))
         OPEN (UNIT=9,FILE=FILTHFILE,STATUS='OLD')
      ENDIF
C 
C FORMAT FOR INFO FILE   LABEL  TYPECH  NUMBER_IN_LIST  BONDEDTO  X  Y  Z
C 
      PRINT *,'MAX. NO. OF ATOMS = ',NATOMS
      DO A=1,NATOMS
        READ (UNIT=9,IOSTAT=IOS,FMT='(A3)') CHECK
        IF (IOS.LT.0) THEN
          PRINT *,'END OF FILE BEFORE ALL INFORMATION SPECIFIED',A
          STOP
        END IF
        IF (CHECK.EQ.'END' .OR. CHECK.EQ.'END' .OR. CHECK.EQ.'END') GOTO 10
        BACKSPACE 9
        READ (UNIT=9,FMT=*) LABEL(A),TYPECH(A),ATNUM(A),BONDEDTO(A),X(A),Y(A),Z(A)
        ZSYM(A)=TYPECH(A)(1:1)
        IF (LABEL(A).EQ.'*') THEN
          CHIRAL(A)=1
C         PRINT *,'FOUND CHIRAL ATOM',A
        END IF
      END DO
10    CONTINUE

      PRINT *,'CO-ORDINATES READ SUCCESSFULLY - BONDEDTO(5)= ',BONDEDTO(5)
      ATOMS=A-1
      PRINT *,'ATOMS= ',ATOMS
C 
C CHECK FOR ANY RINGS AND CLOSE THEM!
C 
      DO A=1,10
       READ (UNIT=9,IOSTAT=IOS,FMT='(A4)') CHECK
       IF (CHECK.NE.'LOOP' .AND. A.EQ.10) GOTO 20
       IF (.NOT.(CHECK.NE.'LOOP')) THEN
          BACKSPACE 9
          READ (UNIT=9,FMT='(A4,7X,I2)') CHECK,RINGS
          GOTO 20
      ENDIF
      END DO
20    CONTINUE

      BEAN=1
      IF (A.LT.10) THEN
       DO I=1,RINGS
C        READ (UNIT=9,FMT='(I3,2X,I3)') LOOPATOM(2*I-1),LOOPATOM(2*I)
        READ (UNIT=9,FMT=*) LOOPATOM(2*I-1),LOOPATOM(2*I)
         B=LOOPATOM(2*I-1)
         C=LOOPATOM(2*I)
         IF (B.EQ.0 .OR. C.EQ.0) THEN
          PRINT *,'NO ATOMS SPECIFIED FOR LOOP BOND'
         END IF
        BONDS(B,C)=1
       BONDS(C,B)=1
       BONDARRAY(BEAN,1)=B
       BONDARRAY(BEAN,2)=C
       BEAN=BEAN+1
       END DO
      ELSE
       PRINT *,'NO LOOP LINE IN CO-ORDINATE FILE'
       STOP
      END IF
C 
C READ CHARGES
C 
      DO A=1,10
       READ (UNIT=9,IOSTAT=IOS,FMT='(A7)') CHECK
       IF (CHECK.NE.'CHARGES' .AND. A.EQ.10) GOTO 30
       IF (.NOT.(CHECK.NE.'CHARGES')) GOTO 30
      END DO
30    CONTINUE

      IF (A.EQ.10) THEN
       PRINT *,'NO CHARGES SPECIFIED'
       STOP
      ELSE 
      DO C=1,ATOMS
C       READ (UNIT=9,IOSTAT=IOS,FMT='(I3,2X,F7.4)') B,VIXEN1
       READ (UNIT=9,IOSTAT=IOS,FMT=*) B,VIXEN1
       PQ(B)=VIXEN1*CRAP
       IF (IOS.LT.0) GOTO 40
      END DO
40    CONTINUE
      END IF
C 
C CONVERTING TYPECH TO TYPE NUMBER
C 
      DO A=1,ATOMS
       CALL TYPENUMBER(TYPECH(A))
       TYPE(A)=ANS
       IF (TYPE(A).EQ.-1) THEN
        PRINT *,'ATOM ',A,'TYPE SPECIFIED INCORRECTLY - ABORTING'
      STOP
       END IF
      END DO
C 
C CONNECTING BACKBONE ATOMS
C 
      DO A=2,ATOMS
       IF (LABEL(A).NE.'B') GOTO 50
       BONDS(A,A-1)=1
       BONDS(A-1,A)=1
      END DO       
50    CONTINUE
C 
C CONNECTING SIDE-CHAIN ATOMS
C 
      DO A=2,ATOMS
       IF (.NOT.(LABEL(A).EQ.'B')) THEN
       BONDS(A,BONDEDTO(A))=1
       BONDS(BONDEDTO(A),A)=1
       BONDARRAY(BEAN,1)=A
       BONDARRAY(BEAN,2)=BONDEDTO(A)
       BEAN=BEAN+1
       ENDIF
      END DO

      BONDNUMBER=BEAN-1

      CALL ONE3

      RETURN
      END

      SUBROUTINE ONE3
      USE MODAMBER
      USE MODAMBER2
       IMPLICIT NONE
C      PRINT *,'SUBROUTINE ONE3'
       DO A=1,ATOMS
         DO B=A+1,ATOMS
           IF (BONDS(A,B) .NE. 1) THEN
             DO C=1,ATOMS
               IF (BONDS(A,C).EQ.1 .AND. BONDS(B,C).EQ.1) THEN
                 ONE_THREE(A,B)=1
                 ONE_THREE(B,A)=1
C                PRINT *,A,B
               END IF
             END DO
           END IF
         END DO
       END DO

      RETURN

      END
      SUBROUTINE APARAMS
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE
      DOUBLE PRECISION VIXEN1, VIXEN2, VIXEN3, VIXEN4
      DOUBLE PRECISION ANGFAC
      PARAMETER (ANGFAC=57.29577951D0)
     
      PRINT *,'READING PARAMETERS'
      HELLCOUNT=0

      OPEN (UNIT=9,FILE='AMBER.DAT',STATUS='OLD')
C 
C READING PARAMETERS
C 
C FORMAT OF PARAMS.DAT FOR BONDS  TYPECHB TYPECHC KR(TYPECHB,TYPECHC) RO(TYPECHB,TYPECHC)
C 
      DO A=1,200
         READ (UNIT=9,IOSTAT=IOS,FMT='(A3)') CHECK
         IF (IOS.LT.0) GOTO 10
         IF (CHECK.EQ.'END' .OR. CHECK.EQ.'END' .OR. CHECK.EQ.'END') GOTO 10
         BACKSPACE 9
         READ (UNIT=9,FMT='(A2,X,A2)') TYPECHB,TYPECHC

         CALL TYPENUMBER(TYPECHB)
         B=ANS
         CALL TYPENUMBER(TYPECHC)
         C=ANS 
         IF (B.EQ.-1 .OR. C.EQ.-1) THEN
          PRINT *,'BOND PARAMETERS ',TYPECHB,'-',TYPECHC,' SPECIFIED INCORRECTLY - ABORTING'
          STOP
         ENDIF

         BACKSPACE 9
         READ (UNIT=9,FMT='(A2,X,A2,2X,F5.0,4X,F6.0)') TYPECHB,TYPECHC,KR(B,C),RO(B,C)
C        WRITE(*,'(3I3,X,A,2X,A,2X,2F15.5)') A,B,C,TYPECHB,TYPECHC,KR(B,C),RO(B,C)
         KR(C,B)=KR(B,C)
         RO(C,B)=RO(B,C)
      ENDDO
10    CONTINUE
C 
C FORMAT OF PARAMS.DAT FOR ANGLES  TYPE1 TYPE2 TYPE3 KT(1,2,3) DEGTO(1,2,3)
C 
      DO A=1,300
         READ (UNIT=9,IOSTAT=IOS,FMT='(A5)') CHECK 
         IF (IOS.LT.0) GOTO 20
         IF (.NOT.(CHECK.EQ.'START' .OR. CHECK.EQ.'START' .OR. CHECK.EQ.'START')) THEN
         IF (.NOT.(CHECK.NE.'START' .AND. CHECK.NE.'START' .AND. CHECK.NE.'START' .AND. A.EQ.1)) THEN
         IF (CHECK.EQ.'END' .OR. CHECK.EQ.'END' .OR. CHECK.EQ.'END') GOTO 20
         BACKSPACE 9
         READ (UNIT=9,FMT='(A2,X,A2,X,A2)') TYPECHB,TYPECHC,TYPECHD

         CALL TYPENUMBER(TYPECHB)
         B=ANS
         CALL TYPENUMBER(TYPECHC)
         C=ANS
         CALL TYPENUMBER(TYPECHD)
         D=ANS

         IF (B.EQ.-1 .OR. C.EQ.-1 .OR. D.EQ.-1) THEN
            PRINT *,'ANGLE ',TYPECHB,'-',TYPECHC,'-',TYPECHD,' SPECIFIED INCORRECTLY - ABORTING'
            STOP
         ENDIF
   
         BACKSPACE 9
         READ (UNIT=9,FMT='(A2,X,A2,X,A2,3X,F5.0,6X,F6.0)') TYPECHB,TYPECHC,TYPECHD,KT(B,C,D),DEGTO(B,C,D)
         TO(B,C,D)=DEGTO(B,C,D)/ANGFAC
         KT(D,C,B)=KT(B,C,D)
         TO(D,C,B)=TO(B,C,D)
         DEGTO(D,C,B)=DEGTO(B,C,D)
         ENDIF
         ENDIF
      ENDDO
20    CONTINUE

      IF (A.EQ.1) THEN
       PRINT *,'ERROR READING ANGLE PARAMETERS - NO START'
       STOP
      ENDIF 
C 
C FORMAT OF PARAMS.DAT FOR TORSIONS  TYPE1 TYPE2 TYPE3 TYPE4 PK IDIVF PN PHASE
C 
      DO A=1,100
         READ (UNIT=9,IOSTAT=IOS,FMT='(A5)') CHECK
         IF (IOS.LT.0) GOTO 30
         IF (.NOT.(CHECK.EQ.'START' .OR. CHECK.EQ.'START' .OR. CHECK.EQ.'START')) THEN
         IF (CHECK.NE.'START' .AND. CHECK.NE.'START' .AND. CHECK.NE.'START' .AND. A.EQ.1) GOTO 30
         IF (CHECK.EQ.'END' .OR. CHECK.EQ.'END' .OR. CHECK.EQ.'END')  GOTO 30
         IF (A.GT.1) BACKSPACE 9
         
         READ (UNIT=9,IOSTAT=IOS,FMT='(A2)') TYPECHB
         IF (TYPECHB.NE.'X ') GOTO 30
         BACKSPACE 9

         READ (UNIT=9,FMT='(A2,X,A2,X,A2,X,A2,3X,F1.0,3X,F6.0,7X,F5.0,13X,F2.0)') TYPECHB,TYPECHC,TYPECHD,
     1        TYPECHE,VIXEN1,VIXEN2,VIXEN3,VIXEN4
         CALL TYPENUMBER(TYPECHC)
           C=ANS
         CALL TYPENUMBER(TYPECHD)
           D=ANS

         GENTORSPARAMS(A-1,1)=C
         GENTORSPARAMS(A-1,2)=D
         GENTORSPARAMS(A-1,3)=VIXEN1 
         GENTORSPARAMS(A-1,4)=VIXEN2
         GENTORSPARAMS(A-1,5)=VIXEN3/ANGFAC
         GENTORSPARAMS(A-1,6)=ABS(VIXEN4)

         ENDIF
      ENDDO
30    CONTINUE
 
      IF (A.EQ.1) THEN
       PRINT *,'ERROR READING TORSION PARAMETERS - NO START'
      STOP
      ENDIF

      BACKSPACE 9
      DO A=1,50
        READ (UNIT=9,IOSTAT=IOS,FMT='(A3)') CHECK
C        PRINT *,"CHECK= ",CHECK
        IF (CHECK.EQ."END") THEN
C        PRINT *,"END OF TORSION PARAMETERS REACHED - A=",A
        GOTO 32
        END IF
        BACKSPACE 9

        READ (UNIT=9,IOSTAT=IOS,FMT='(A2,X,A2,X,A2,X,A2,3X,F1.0,3X,F6.0,7X,F5.0,12X,F3.0)') TYPECHB,TYPECHC,
     1       TYPECHD,TYPECHE,VIXEN1,VIXEN2,VIXEN3,VIXEN4
        CALL TYPENUMBER(TYPECHB)
          B=ANS
        CALL TYPENUMBER(TYPECHC)
          C=ANS
        CALL TYPENUMBER(TYPECHD)
          D=ANS
        CALL TYPENUMBER(TYPECHE)
          E=ANS

        SPECTORSPARAMS(A,1)=B
        SPECTORSPARAMS(A,2)=C
        SPECTORSPARAMS(A,3)=D
        SPECTORSPARAMS(A,4)=E
        SPECTORSPARAMS(A,5)=VIXEN1
        SPECTORSPARAMS(A,6)=VIXEN2
        SPECTORSPARAMS(A,7)=VIXEN3/ANGFAC
        SPECTORSPARAMS(A,8)=ABS(VIXEN4)

        IF (VIXEN4.GT.0) THEN
         SPECTORSPARAMS(A,9)=0.0
         SPECTORSPARAMS(A,10)=0.0
         SPECTORSPARAMS(A,11)=0.0
         SPECTORSPARAMS(A,12)=0.0
         SPECTORSPARAMS(A,13)=0.0
         SPECTORSPARAMS(A,14)=0.0
         GOTO 31
        END IF

        READ (UNIT=9,IOSTAT=IOS,FMT='(A2,X,A2,X,A2,X,A2,3X,F1.0,3X,F6.0,7X,F5.0,12X,F2.0)') TYPECHB,TYPECHC,
     1      TYPECHD,TYPECHE,VIXEN1,VIXEN2,VIXEN3,VIXEN4

        SPECTORSPARAMS(A,9)=VIXEN2
        SPECTORSPARAMS(A,10)=VIXEN3/ANGFAC
        SPECTORSPARAMS(A,11)=ABS(VIXEN4)
       
        IF (VIXEN4.GT.0) THEN
          SPECTORSPARAMS(A,12)=0.0
          SPECTORSPARAMS(A,13)=0.0
          SPECTORSPARAMS(A,14)=0.0
          GOTO 31
        END IF

        READ (UNIT=9,IOSTAT=IOS,FMT='(A2,X,A2,X,A2,X,A2,3X,F1.0,3X,F6.0,7X,F5.0,12X,F2.0)') TYPECHB,TYPECHC,
     1    TYPECHD,TYPECHE,VIXEN1,VIXEN2,VIXEN3,VIXEN4

        SPECTORSPARAMS(A,12)=VIXEN2
        SPECTORSPARAMS(A,13)=VIXEN3/ANGFAC
        SPECTORSPARAMS(A,14)=VIXEN4

31     CONTINUE
       END DO
32     CONTINUE
C 
C FORMAT FOR IMPROPER TORSION PARAMETERS    TYPECHA  TYPECHB  TYPECHC  TYPECHD  IPK  IPHASE  IPN
C 
      DO A=1,100
       READ (UNIT=9,IOSTAT=IOS,FMT='(A9)') CHECK
       IF (IOS.LT.0) THEN
        PRINT *,'END OF PARAMETERS FILE WITHOUT IMPROPER TORSION PARAMETERS'
        STOP
       ENDIF
       IF (.NOT.(CHECK.NE.'IMPROPERS' .AND. CHECK.NE.'IMPROPERS' .AND. CHECK.NE.'IMPROPERS')) GOTO 40
      ENDDO
40    CONTINUE

      IF (A.GE.100) THEN 
       PRINT *,'NO IMPROPER TORSION PARAMETERS SPECIFIED'
      STOP
      ENDIF

      BEAN=0
      DO A=1,100
       READ (UNIT=9,IOSTAT=IOS,FMT='(A2,X,A2,X,A2,X,A2)') TYPECHB,TYPECHC,TYPECHD,TYPECHE
       IF (IOS.LT.0) THEN
        PRINT *,'END OF FILE DURING IMPROPER TORSION PARAMETERS'
      STOP
       ENDIF
      IF (TYPECHB.EQ.'X ' .AND. TYPECHC.EQ.'X ') THEN
        BACKSPACE 9 
        READ (UNIT=9,FMT='(A2,X,A2,X,A2,X,A2)') TYPECHB,TYPECHC,TYPECHD,TYPECHE
        CALL TYPENUMBER(TYPECHD)
          D=ANS
        CALL TYPENUMBER(TYPECHE)
          E=ANS
 
        GENIMPPARAMS(A-BEAN,1)=D
        GENIMPPARAMS(A-BEAN,2)=E

        BACKSPACE 9
        READ (UNIT=9,FMT='(A2,X,A2,X,A2,X,A2,9X,F4.0,9X,F4.0,10X,F1.0)') TYPECHB,TYPECHC,TYPECHD,TYPECHE,
     1       VIXEN1,VIXEN2,VIXEN3
        GENIMPPARAMS(A-BEAN,3)=VIXEN1
        GENIMPPARAMS(A-BEAN,4)=VIXEN2/ANGFAC
        GENIMPPARAMS(A-BEAN,5)=VIXEN3

       ELSE IF (TYPECHB.EQ.'X ' .AND. TYPECHC.NE.'X ') THEN
        BEAN=BEAN+1
        BACKSPACE 9
        READ (UNIT=9,FMT='(A2,X,A2,X,A2,X,A2)') TYPECHB,TYPECHC,TYPECHD,TYPECHE
         CALL TYPENUMBER(TYPECHC)
           C=ANS
         CALL TYPENUMBER(TYPECHD)
           D=ANS
         CALL TYPENUMBER(TYPECHE)
           E=ANS

        BACKSPACE 9
        READ (UNIT=9,FMT='(A2,X,A2,X,A2,X,A2,9X,F4.0,9X,F4.0,10X,F1.0)') TYPECHB,TYPECHC,TYPECHD,TYPECHE,
     1       VIXEN1,VIXEN2,VIXEN3

        MIDIMPPARAMS(BEAN,1)=C
        MIDIMPPARAMS(BEAN,2)=D        
        MIDIMPPARAMS(BEAN,3)=E
        MIDIMPPARAMS(BEAN,4)=VIXEN1
        MIDIMPPARAMS(BEAN,5)=VIXEN2/ANGFAC
        MIDIMPPARAMS(BEAN,6)=VIXEN3

       ELSE IF (TYPECHB.EQ.'  ') THEN
          GOTO 50

       ELSE IF (TYPECHB.EQ.'EN' .OR. TYPECHB.EQ.'EN' .OR. TYPECHB.EQ.'EN') THEN
          GOTO 60
       ELSE 
          GOTO 60

       ENDIF
50     CONTINUE
      ENDDO   
60    CONTINUE

      BACKSPACE 9
      DO A=1,20
        READ (UNIT=9,FMT='(A2,X,A2,X,A2,X,A2,9X,F4.0,9X,F4.0,10X,F1.0)') TYPECHB,TYPECHC,TYPECHD,TYPECHE,
     1       VIXEN1,VIXEN2,VIXEN3
        IF (TYPECHB.EQ.'EN') GOTO 65
        CALL TYPENUMBER(TYPECHB)
          B=ANS
        CALL TYPENUMBER(TYPECHC)
          C=ANS
        CALL TYPENUMBER(TYPECHD)
          D=ANS
        CALL TYPENUMBER(TYPECHE)
          E=ANS

        SPECIMPPARAMS(A,1)=B
        SPECIMPPARAMS(A,2)=C
        SPECIMPPARAMS(A,3)=D
        SPECIMPPARAMS(A,4)=E
        SPECIMPPARAMS(A,5)=VIXEN1
        SPECIMPPARAMS(A,6)=VIXEN2/ANGFAC
        SPECIMPPARAMS(A,7)=VIXEN3

      END DO
65    CONTINUE

      DO A=1,100
       READ (UNIT=9,IOSTAT=IOS,FMT='(A3)') CHECK
       IF (IOS.LT.0) THEN
        PRINT *,'END OF PARAMETER FILE WITHOUT VDW PARAMETERS'
        STOP
       ENDIF
       IF (.NOT.(CHECK.NE.'VDW' .AND. CHECK.NE.'VDW')) GOTO 70
      ENDDO
70    CONTINUE
  
      IF (A.LT.100) THEN
       DO A=1,42
        READ (UNIT=9,IOSTAT=IOS,FMT='(A2)') TYPECHB
        IF (TYPECHB.EQ.'EN') GOTO 80
        IF (IOS.LT.0) GOTO 80
        CALL TYPENUMBER(TYPECHB)
        B=ANS
        BACKSPACE 9
        READ (UNIT=9,FMT='(A2,10X,F6.4,2X,F6.4)') TYPECHB,VDWR(B),VDWE(B)
       ENDDO
80     CONTINUE
      ELSE
       PRINT *,'UNABLE TO FIND VDW PARAMETERS'
      ENDIF

      PRINT *,'PARAMETERS READ'
    
      RETURN
      END

      SUBROUTINE TYPENUMBER(TYPECHAR)
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE
      CHARACTER(LEN=2)  TYPECHAR

      IF (TYPECHAR.EQ.'C ') THEN
         ANS=1
      ELSE IF (TYPECHAR.EQ.'CA') THEN
         ANS=2
      ELSE IF (TYPECHAR.EQ.'CB') THEN
         ANS=3
      ELSE IF (TYPECHAR.EQ.'CC') THEN
         ANS=4
      ELSE IF (TYPECHAR.EQ.'CK') THEN
         ANS=5
      ELSE IF (TYPECHAR.EQ.'CM') THEN
         ANS=6
      ELSE IF (TYPECHAR.EQ.'CN') THEN
         ANS=7
      ELSE IF (TYPECHAR.EQ.'CQ') THEN
         ANS=8
      ELSE IF (TYPECHAR.EQ.'CR') THEN
         ANS=9
      ELSE IF (TYPECHAR.EQ.'CT') THEN
         ANS=10
      ELSE IF (TYPECHAR.EQ.'CV') THEN
         ANS=11
      ELSE IF (TYPECHAR.EQ.'CW') THEN
         ANS=12
      ELSE IF (TYPECHAR.EQ.'C*') THEN
         ANS=13
      ELSE IF (TYPECHAR.EQ.'C0') THEN
         ANS=14
      ELSE IF (TYPECHAR.EQ.'F ') THEN
         ANS=15
      ELSE IF (TYPECHAR.EQ.'H ') THEN
         ANS=16
      ELSE IF (TYPECHAR.EQ.'HC') THEN
         ANS=17
      ELSE IF (TYPECHAR.EQ.'H1') THEN
         ANS=18
      ELSE IF (TYPECHAR.EQ.'H2') THEN
         ANS=19
      ELSE IF (TYPECHAR.EQ.'H3') THEN
         ANS=20
      ELSE IF (TYPECHAR.EQ.'HA') THEN
         ANS=21
      ELSE IF (TYPECHAR.EQ.'H4') THEN
         ANS=22
      ELSE IF (TYPECHAR.EQ.'H5') THEN
         ANS=23
      ELSE IF (TYPECHAR.EQ.'HO') THEN
         ANS=24
      ELSE IF (TYPECHAR.EQ.'HS') THEN
         ANS=25
      ELSE IF (TYPECHAR.EQ.'HW') THEN
         ANS=26
      ELSE IF (TYPECHAR.EQ.'HP') THEN
         ANS=27
      ELSE IF (TYPECHAR.EQ.'N ') THEN
         ANS=28
      ELSE IF (TYPECHAR.EQ.'NA') THEN
         ANS=29
      ELSE IF (TYPECHAR.EQ.'NB') THEN
         ANS=30
      ELSE IF (TYPECHAR.EQ.'NC') THEN
         ANS=31
      ELSE IF (TYPECHAR.EQ.'N2') THEN
         ANS=32
      ELSE IF (TYPECHAR.EQ.'N3') THEN
         ANS=33
      ELSE IF (TYPECHAR.EQ.'N*') THEN
         ANS=34
      ELSE IF (TYPECHAR.EQ.'O ') THEN
         ANS=35
      ELSE IF (TYPECHAR.EQ.'OW') THEN
         ANS=36
      ELSE IF (TYPECHAR.EQ.'OH') THEN
         ANS=37
      ELSE IF (TYPECHAR.EQ.'OS') THEN
         ANS=38
      ELSE IF (TYPECHAR.EQ.'O2') THEN
         ANS=39
      ELSE IF (TYPECHAR.EQ.'P ') THEN
         ANS=40
      ELSE IF (TYPECHAR.EQ.'S ') THEN
         ANS=41
      ELSE IF (TYPECHAR.EQ.'SH') THEN
         ANS=42
      ELSE 
         ANS=-1
      ENDIF

      RETURN
      END
      SUBROUTINE SECONDDERIVS
      USE COMMONS
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE
      DOUBLE PRECISION      VIXEN1,VIXEN2,VIXEN3,VIXEN4,VIXEN5,VIXEN6,VIXEN7,VIXEN8,VIXEN9,JKL
      DOUBLE PRECISION      TORSCRAP(3*NATOMS)
      DOUBLE PRECISION      MUCRAP, LAMBDACRAP, RBC, RBC2, RBC4, RBC6, XE, YE, ZE, RAB, RAB5, RAB2, RAB3, RAB6
      DOUBLE PRECISION      KBOND, RBOND, RAB4, ST, CT, XAB2, YAB2, ZAB2, XAB, YAB, ZAB

      PRINT *,'IN SECONDDERIVS' 
      HELLCOUNT= HELLCOUNT+1
C
C INITIALISE!!
C
      DO A=1,3*ATOMS
         TORSCRAP(A)=0.0D0
         DO B=1,3*ATOMS
            BONDHELL(A,B)=0.0D0
            ANGLEHELL(A,B)=0.0D0
            TORSHELL(A,B)=0.0D0
            IMPHELL(A,B)=0.0D0
            QHELL(A,B)=0.0D0
            VDWHELL(A,B)=0.0D0
         END DO
      END DO

      DO A=1,ATOMS
         DO B=A+1,ATOMS
            RAB2=R(A,B)**2
            RAB3=R(A,B)**3
            RAB4=RAB2**2
            RAB5=RAB2*RAB3
            RAB6=RAB3**2
            XAB=X(A)-X(B)
            YAB=Y(A)-Y(B)
            ZAB=Z(A)-Z(B)
            XAB2=XAB**2
            YAB2=YAB**2
            ZAB2=ZAB**2

            IF (BONDS(A,B).NE.1) GOTO 201

            KBOND=KR(TYPE(A),TYPE(B))
            RBOND=RO(TYPE(A),TYPE(B))
            VIXEN1=RBOND/R(A,B)
            VIXEN2=((X(A)-X(B))**2)*RBOND/(RAB3)
            JKL= 2.0D0*KBOND*(1.0D0-VIXEN1+VIXEN2)
C
C DXADXA,DXBDXB,DXADXB,DXBDXA
C
            BONDHELL(3*A-2,3*A-2)= BONDHELL(3*A-2,3*A-2)+ JKL
            BONDHELL(3*B-2,3*B-2)= BONDHELL(3*B-2,3*B-2)+ JKL
            BONDHELL(3*A-2,3*B-2)= BONDHELL(3*A-2,3*B-2)- JKL 
            BONDHELL(3*B-2,3*A-2)= BONDHELL(3*B-2,3*A-2)- JKL 

            VIXEN2=((Y(A)-Y(B))**2)*RBOND/(RAB3)
            JKL= 2.0D0*KBOND*(1.0D0-VIXEN1+VIXEN2)
C
C DYADYA,DYBDYB,DYADYB,DYBDYA
C
            BONDHELL(3*A-1,3*A-1)= BONDHELL(3*A-1,3*A-1)+ JKL
            BONDHELL(3*B-1,3*B-1)= BONDHELL(3*B-1,3*B-1)+ JKL
            BONDHELL(3*A-1,3*B-1)= BONDHELL(3*A-1,3*B-1)- JKL
            BONDHELL(3*B-1,3*A-1)= BONDHELL(3*B-1,3*A-1)- JKL

            VIXEN2=((Z(A)-Z(B))**2)*RBOND/(RAB3)
            JKL= 2.0D0*KBOND*(1.0D0-VIXEN1+VIXEN2)
C
C DZADZA,DZBDZB,DZADZB,DZBDZA
C
            BONDHELL(3*A,3*A)= BONDHELL(3*A,3*A)+ JKL
            BONDHELL(3*B,3*B)= BONDHELL(3*B,3*B)+ JKL
            BONDHELL(3*A,3*B)= BONDHELL(3*A,3*B)- JKL
            BONDHELL(3*B,3*A)= BONDHELL(3*B,3*A)- JKL

            VIXEN1=2.0D0*KBOND*(X(A)-X(B))*(Y(A)-Y(B))*RBOND/RAB3
C
C DXADYA,DYADXA,DXBDYB,DYBDXB
C
            BONDHELL(3*A-2,3*A-1)= BONDHELL(3*A-2,3*A-1)+ VIXEN1
            BONDHELL(3*A-1,3*A-2)= BONDHELL(3*A-1,3*A-2)+ VIXEN1
            BONDHELL(3*B-2,3*B-1)= BONDHELL(3*B-2,3*B-1)+ VIXEN1
            BONDHELL(3*B-1,3*B-2)= BONDHELL(3*B-1,3*B-2)+ VIXEN1

            VIXEN1=2.0D0*KBOND*(X(A)-X(B))*(Z(A)-Z(B))*RBOND/RAB3
C
C DXADZA,DZADXA,DXBDZB,DZBDXB
C
            BONDHELL(3*A-2,3*A)= BONDHELL(3*A-2,3*A)+ VIXEN1
            BONDHELL(3*A,3*A-2)= BONDHELL(3*A,3*A-2)+ VIXEN1
            BONDHELL(3*B-2,3*B)= BONDHELL(3*B-2,3*B)+ VIXEN1
            BONDHELL(3*B,3*B-2)= BONDHELL(3*B,3*B-2)+ VIXEN1

            VIXEN1=2.0D0*KBOND*(Z(A)-Z(B))*(Y(A)-Y(B))*RBOND/RAB3
C
C DYADZA,DZADYA,DYBDZB,DZBDYB
C
            BONDHELL(3*A-1,3*A)= BONDHELL(3*A-1,3*A)+ VIXEN1
            BONDHELL(3*A,3*A-1)= BONDHELL(3*A,3*A-1)+ VIXEN1
            BONDHELL(3*B-1,3*B)= BONDHELL(3*B-1,3*B)+ VIXEN1
            BONDHELL(3*B,3*B-1)= BONDHELL(3*B,3*B-1)+ VIXEN1

            VIXEN1=-2.0D0*KBOND*(X(A)-X(B))*(Y(A)-Y(B))*RBOND/RAB3
C
C DXADYB,DYBDXA,DXBDYA,DYADXB
C
            BONDHELL(3*A-2,3*B-1)= BONDHELL(3*A-2,3*B-1)+ VIXEN1
            BONDHELL(3*B-1,3*A-2)= BONDHELL(3*B-1,3*A-2)+ VIXEN1
            BONDHELL(3*B-2,3*A-1)= BONDHELL(3*B-2,3*A-1)+ VIXEN1
            BONDHELL(3*A-1,3*B-2)= BONDHELL(3*A-1,3*B-2)+ VIXEN1

            VIXEN1=-2.0D0*KBOND*(X(A)-X(B))*(Z(A)-Z(B))*RBOND/RAB3
C
C DXADZB,DZBDXA,DXBDZA,DZADXB
C
            BONDHELL(3*A-2,3*B)= BONDHELL(3*A-2,3*B)+ VIXEN1
            BONDHELL(3*B,3*A-2)= BONDHELL(3*B,3*A-2)+ VIXEN1
            BONDHELL(3*B-2,3*A)= BONDHELL(3*B-2,3*A)+ VIXEN1
            BONDHELL(3*A,3*B-2)= BONDHELL(3*A,3*B-2)+ VIXEN1

            VIXEN1=-2.0D0*KBOND*(Y(A)-Y(B))*(Z(A)-Z(B))*RBOND/RAB3
C
C DYADZB,DZBDYA,DYBDZA,DZADYB
C
            BONDHELL(3*A-1,3*B)= BONDHELL(3*A-1,3*B)+ VIXEN1
            BONDHELL(3*B,3*A-1)= BONDHELL(3*B,3*A-1)+ VIXEN1
            BONDHELL(3*B-1,3*A)= BONDHELL(3*B-1,3*A)+ VIXEN1
            BONDHELL(3*A,3*B-1)= BONDHELL(3*A,3*B-1)+ VIXEN1

201      CONTINUE
C
C ELECTROSTATIC DERIVATIVES
C
            IF (FAKEWATER) THEN
               VIXEN2=(1.0D0-DBLE(BONDS(A,B)+ONE_THREE(A,B)))/(1.0D0+0.2D0*ONE_FOUR(A,B))
C
C DXADXA,DXBDXB,DXADXB,DXBDXA
C
               QHELL(3*A-2,3*A-2)= QHELL(3*A-2,3*A-2)+ VIXEN2*2.0D0*PQ(A)*PQ(B)*(4.0D0*XAB2/RAB2-1.0D0)/RAB4
               QHELL(3*B-2,3*B-2)= QHELL(3*B-2,3*B-2)+ VIXEN2*2.0D0*PQ(A)*PQ(B)*(4.0D0*XAB2/RAB2-1.0D0)/RAB4
               QHELL(3*A-2,3*B-2)= QHELL(3*A-2,3*B-2)- VIXEN2*2.0D0*PQ(A)*PQ(B)*(4.0D0*XAB2/RAB2-1.0D0)/RAB4
               QHELL(3*B-2,3*A-2)= QHELL(3*B-2,3*A-2)- VIXEN2*2.0D0*PQ(A)*PQ(B)*(4.0D0*XAB2/RAB2-1.0D0)/RAB4
C
C DYADYA,DYBDYB,DYADYB,DYBDYA
C
               QHELL(3*A-1,3*A-1)= QHELL(3*A-1,3*A-1)+ VIXEN2*2.0D0*PQ(A)*PQ(B)*(4.0D0*YAB2/RAB2-1.0D0)/RAB4
               QHELL(3*B-1,3*B-1)= QHELL(3*B-1,3*B-1)+ VIXEN2*2.0D0*PQ(A)*PQ(B)*(4.0D0*YAB2/RAB2-1.0D0)/RAB4
               QHELL(3*A-1,3*B-1)= QHELL(3*A-1,3*B-1)- VIXEN2*2.0D0*PQ(A)*PQ(B)*(4.0D0*YAB2/RAB2-1.0D0)/RAB4
               QHELL(3*B-1,3*A-1)= QHELL(3*B-1,3*A-1)- VIXEN2*2.0D0*PQ(A)*PQ(B)*(4.0D0*YAB2/RAB2-1.0D0)/RAB4
C
C DZADZA,DZBDZB,DZADZB,DZBDZA
C
               QHELL(3*A,3*A)= QHELL(3*A,3*A)+ VIXEN2*2.0D0*PQ(A)*PQ(B)*(4.0D0*ZAB2/RAB2-1.0D0)/RAB4
               QHELL(3*B,3*B)= QHELL(3*B,3*B)+ VIXEN2*2.0D0*PQ(A)*PQ(B)*(4.0D0*ZAB2/RAB2-1.0D0)/RAB4
               QHELL(3*A,3*B)= QHELL(3*A,3*B)- VIXEN2*2.0D0*PQ(A)*PQ(B)*(4.0D0*ZAB2/RAB2-1.0D0)/RAB4
               QHELL(3*B,3*A)= QHELL(3*B,3*A)- VIXEN2*2.0D0*PQ(A)*PQ(B)*(4.0D0*ZAB2/RAB2-1.0D0)/RAB4
C
C DXADYA,DYADXA,DXBDYB,DYBDXB
C
               QHELL(3*A-2,3*A-1)= QHELL(3*A-2,3*A-1)+ VIXEN2*8.0D0*PQ(A)*PQ(B)*XAB*YAB/RAB6
               QHELL(3*A-1,3*A-2)= QHELL(3*A-1,3*A-2)+ VIXEN2*8.0D0*PQ(A)*PQ(B)*XAB*YAB/RAB6
               QHELL(3*B-2,3*B-1)= QHELL(3*B-2,3*B-1)+ VIXEN2*8.0D0*PQ(A)*PQ(B)*XAB*YAB/RAB6
               QHELL(3*B-1,3*B-2)= QHELL(3*B-1,3*B-2)+ VIXEN2*8.0D0*PQ(A)*PQ(B)*XAB*YAB/RAB6
C
C DXADYB,DYBDXA,DXBDYA,DYADXB
C
               QHELL(3*A-2,3*B-1)= QHELL(3*A-2,3*B-1)- VIXEN2*8.0D0*PQ(A)*PQ(B)*XAB*YAB/RAB6
               QHELL(3*B-1,3*A-2)= QHELL(3*B-1,3*A-2)- VIXEN2*8.0D0*PQ(A)*PQ(B)*XAB*YAB/RAB6
               QHELL(3*B-2,3*A-1)= QHELL(3*B-2,3*A-1)- VIXEN2*8.0D0*PQ(A)*PQ(B)*XAB*YAB/RAB6
               QHELL(3*A-1,3*B-2)= QHELL(3*A-1,3*B-2)- VIXEN2*8.0D0*PQ(A)*PQ(B)*XAB*YAB/RAB6
C
C DXADZA,DZADXA,DXBDZB,DZBDXB
C
               QHELL(3*A-2,3*A)= QHELL(3*A-2,3*A)+ VIXEN2*8.0D0*PQ(A)*PQ(B)*XAB*ZAB/RAB6
               QHELL(3*A,3*A-2)= QHELL(3*A,3*A-2)+ VIXEN2*8.0D0*PQ(A)*PQ(B)*XAB*ZAB/RAB6
               QHELL(3*B-2,3*B)= QHELL(3*B-2,3*B)+ VIXEN2*8.0D0*PQ(A)*PQ(B)*XAB*ZAB/RAB6
               QHELL(3*B,3*B-2)= QHELL(3*B,3*B-2)+ VIXEN2*8.0D0*PQ(A)*PQ(B)*XAB*ZAB/RAB6
C
C DXADZB,DZBDXA,DXBDZA,DZADXB
C
               QHELL(3*A-2,3*B)= QHELL(3*A-2,3*B)- VIXEN2*8.0D0*PQ(A)*PQ(B)*XAB*ZAB/RAB6
               QHELL(3*B,3*A-2)= QHELL(3*B,3*A-2)- VIXEN2*8.0D0*PQ(A)*PQ(B)*XAB*ZAB/RAB6
               QHELL(3*B-2,3*A)= QHELL(3*B-2,3*A)- VIXEN2*8.0D0*PQ(A)*PQ(B)*XAB*ZAB/RAB6
               QHELL(3*A,3*B-2)= QHELL(3*A,3*B-2)- VIXEN2*8.0D0*PQ(A)*PQ(B)*XAB*ZAB/RAB6
C
C DZADYA,DYADZA,DZBDYB,DYBDZB
C
               QHELL(3*A,3*A-1)= QHELL(3*A,3*A-1)+ VIXEN2*8.0D0*PQ(A)*PQ(B)*ZAB*YAB/RAB6
               QHELL(3*A-1,3*A)= QHELL(3*A-1,3*A)+ VIXEN2*8.0D0*PQ(A)*PQ(B)*ZAB*YAB/RAB6
               QHELL(3*B,3*B-1)= QHELL(3*B,3*B-1)+ VIXEN2*8.0D0*PQ(A)*PQ(B)*ZAB*YAB/RAB6
               QHELL(3*B-1,3*B)= QHELL(3*B-1,3*B)+ VIXEN2*8.0D0*PQ(A)*PQ(B)*ZAB*YAB/RAB6
C
C DZADYB,DYBDZA,DZBDYA,DYADZB
C
               QHELL(3*A,3*B-1)= QHELL(3*A,3*B-1)- VIXEN2*8.0D0*PQ(A)*PQ(B)*ZAB*YAB/RAB6
               QHELL(3*B-1,3*A)= QHELL(3*B-1,3*A)- VIXEN2*8.0D0*PQ(A)*PQ(B)*ZAB*YAB/RAB6
               QHELL(3*B,3*A-1)= QHELL(3*B,3*A-1)- VIXEN2*8.0D0*PQ(A)*PQ(B)*ZAB*YAB/RAB6
               QHELL(3*A-1,3*B)= QHELL(3*A-1,3*B)- VIXEN2*8.0D0*PQ(A)*PQ(B)*ZAB*YAB/RAB6


            ELSE
               VIXEN2= (1.0D0-DBLE(BONDS(A,B)+ONE_THREE(A,B)))*PQ(A)*PQ(B)/
     1                 ((1.0D0+0.2D0*ONE_FOUR(A,B))*DIELEC*RAB5)
C
C DXADXA,DXBDXB,DXADXB,DXBDXA
C
               VIXEN1= (3.0D0*XAB**2 - RAB2)*VIXEN2
               QHELL(3*A-2,3*A-2)= QHELL(3*A-2,3*A-2)+ VIXEN1
               QHELL(3*B-2,3*B-2)= QHELL(3*B-2,3*B-2)+ VIXEN1
               QHELL(3*A-2,3*B-2)= QHELL(3*A-2,3*B-2)- VIXEN1
               QHELL(3*B-2,3*A-2)= QHELL(3*B-2,3*A-2)- VIXEN1
C
C DYADYA,DYBDYB,DYADYB,DYBDYA
C
               VIXEN1= (3.0D0*YAB**2 - RAB2)*VIXEN2
               QHELL(3*A-1,3*A-1)= QHELL(3*A-1,3*A-1)+ VIXEN1
               QHELL(3*B-1,3*B-1)= QHELL(3*B-1,3*B-1)+ VIXEN1
               QHELL(3*A-1,3*B-1)= QHELL(3*A-1,3*B-1)- VIXEN1
               QHELL(3*B-1,3*A-1)= QHELL(3*B-1,3*A-1)- VIXEN1
C
C DZADZA,DZBDZB,DZADZB,DZBDZA
C
               VIXEN1= (3.0D0*ZAB**2 - RAB2)*VIXEN2
               QHELL(3*A,3*A)= QHELL(3*A,3*A)+ VIXEN1
               QHELL(3*B,3*B)= QHELL(3*B,3*B)+ VIXEN1
               QHELL(3*A,3*B)= QHELL(3*A,3*B)- VIXEN1
               QHELL(3*B,3*A)= QHELL(3*B,3*A)- VIXEN1

               VIXEN7= VIXEN2*3.0D0*XAB*YAB
               VIXEN8= VIXEN2*3.0D0*XAB*ZAB
               VIXEN9= VIXEN2*3.0D0*YAB*ZAB
C
C DXADYA,DYADXA,DXBDYB,DYBDXB
C
               QHELL(3*A-2,3*A-1)= QHELL(3*A-2,3*A-1)+ VIXEN7
               QHELL(3*A-1,3*A-2)= QHELL(3*A-1,3*A-2)+ VIXEN7
               QHELL(3*B-2,3*B-1)= QHELL(3*B-2,3*B-1)+ VIXEN7
               QHELL(3*B-1,3*B-2)= QHELL(3*B-1,3*B-2)+ VIXEN7
C
C DXADZA,DZADXA,DXBDZB,DZBDXB
C
               QHELL(3*A-2,3*A)= QHELL(3*A-2,3*A)+ VIXEN8
               QHELL(3*A,3*A-2)= QHELL(3*A,3*A-2)+ VIXEN8
               QHELL(3*B-2,3*B)= QHELL(3*B-2,3*B)+ VIXEN8
               QHELL(3*B,3*B-2)= QHELL(3*B,3*B-2)+ VIXEN8
C
C DYADZA,DZADYA,DYBDZB,DZBDYB
C
               QHELL(3*A-1,3*A)= QHELL(3*A-1,3*A)+ VIXEN9
               QHELL(3*A,3*A-1)= QHELL(3*A,3*A-1)+ VIXEN9
               QHELL(3*B-1,3*B)= QHELL(3*B-1,3*B)+ VIXEN9
               QHELL(3*B,3*B-1)= QHELL(3*B,3*B-1)+ VIXEN9
C
C DXADYB,DYBDXA,DXBDYA,DYADXB
C
               QHELL(3*A-2,3*B-1)= QHELL(3*A-2,3*B-1)- VIXEN7
               QHELL(3*B-1,3*A-2)= QHELL(3*B-1,3*A-2)- VIXEN7
               QHELL(3*B-2,3*A-1)= QHELL(3*B-2,3*A-1)- VIXEN7
               QHELL(3*A-1,3*B-2)= QHELL(3*A-1,3*B-2)- VIXEN7
C
C DXADZB,DZBDXA,DXBDZA,DZADXB
C
               QHELL(3*A-2,3*B)= QHELL(3*A-2,3*B)- VIXEN8
               QHELL(3*B,3*A-2)= QHELL(3*B,3*A-2)- VIXEN8
               QHELL(3*B-2,3*A)= QHELL(3*B-2,3*A)- VIXEN8
               QHELL(3*A,3*B-2)= QHELL(3*A,3*B-2)- VIXEN8
C
C DYADZB,DZBDYA,DYBDZA,DZADYB
C            
               QHELL(3*A-1,3*B)= QHELL(3*A-1,3*B)- VIXEN9
               QHELL(3*B,3*A-1)= QHELL(3*B,3*A-1)- VIXEN9
               QHELL(3*B-1,3*A)= QHELL(3*B-1,3*A)- VIXEN9
               QHELL(3*A,3*B-1)= QHELL(3*A,3*B-1)- VIXEN9
            END IF

C
C VDW DERIVATIVES
C
            VIXEN2= (1.0D0-DBLE(BONDS(A,B)+ONE_THREE(A,B)))/(1.0D0+DBLE(ONE_FOUR(A,B)))
C
C DXADXA,DXBDXB,DXADXB,DXBDXA
C
            VIXEN1=24.0D0*(7.0D0*(VDWA(A,B)/RAB4**4)-2.0D0*(VDWB(A,B)/(RAB4**2*RAB2)))
            VIXEN3=6.0D0*(VDWB(A,B)/RAB4**2-2.0D0*(VDWA(A,B)/(RAB4**3*RAB2)))
            JKL= VIXEN2*(VIXEN1*XAB2+VIXEN3)
            VDWHELL(3*A-2,3*A-2)= VDWHELL(3*A-2,3*A-2)+ JKL
            VDWHELL(3*B-2,3*B-2)= VDWHELL(3*B-2,3*B-2)+ JKL
            VDWHELL(3*A-2,3*B-2)= VDWHELL(3*A-2,3*B-2)- JKL
            VDWHELL(3*B-2,3*A-2)= VDWHELL(3*B-2,3*A-2)- JKL
C
C DYADYA,DYBDYB,DYADYB,DYBDYA
C
            JKL= VIXEN2*(VIXEN1*YAB2+VIXEN3)
            VDWHELL(3*A-1,3*A-1)= VDWHELL(3*A-1,3*A-1)+ JKL
            VDWHELL(3*B-1,3*B-1)= VDWHELL(3*B-1,3*B-1)+ JKL
            VDWHELL(3*A-1,3*B-1)= VDWHELL(3*A-1,3*B-1)- JKL
            VDWHELL(3*B-1,3*A-1)= VDWHELL(3*B-1,3*A-1)- JKL
C
C DZADZA,DZBDZB,DZADZB,DZBDZA
C
            JKL= VIXEN2*(VIXEN1*ZAB2+VIXEN3)
            VDWHELL(3*A,3*A)= VDWHELL(3*A,3*A)+ JKL
            VDWHELL(3*B,3*B)= VDWHELL(3*B,3*B)+ JKL
            VDWHELL(3*A,3*B)= VDWHELL(3*A,3*B)- JKL
            VDWHELL(3*B,3*A)= VDWHELL(3*B,3*A)- JKL
C
C DXADYA,DYADXA,DXBDYB,DYBDXB
C
            JKL= VIXEN2*VIXEN1*YAB*XAB
            VDWHELL(3*A-2,3*A-1)= VDWHELL(3*A-2,3*A-1)+ JKL
            VDWHELL(3*A-1,3*A-2)= VDWHELL(3*A-1,3*A-2)+ JKL
            VDWHELL(3*B-2,3*B-1)= VDWHELL(3*B-2,3*B-1)+ JKL
            VDWHELL(3*B-1,3*B-2)= VDWHELL(3*B-1,3*B-2)+ JKL
C
C DXADYB,DYBDXA,DXBDYA,DYADXB
C
            VDWHELL(3*A-2,3*B-1)= VDWHELL(3*A-2,3*B-1)- JKL
            VDWHELL(3*B-1,3*A-2)= VDWHELL(3*B-1,3*A-2)- JKL
            VDWHELL(3*B-2,3*A-1)= VDWHELL(3*B-2,3*A-1)- JKL
            VDWHELL(3*A-1,3*B-2)= VDWHELL(3*A-1,3*B-2)- JKL
C
C DXADZA,DZADXA,DXBDZB,DZBDXB
C
            JKL= VIXEN2*VIXEN1*XAB*ZAB
            VDWHELL(3*A-2,3*A)= VDWHELL(3*A-2,3*A)+ JKL
            VDWHELL(3*A,3*A-2)= VDWHELL(3*A,3*A-2)+ JKL
            VDWHELL(3*B-2,3*B)= VDWHELL(3*B-2,3*B)+ JKL
            VDWHELL(3*B,3*B-2)= VDWHELL(3*B,3*B-2)+ JKL
C
C DXADZB,DZBDXA,DXBDZA,DZADXB
C
            VDWHELL(3*A-2,3*B)= VDWHELL(3*A-2,3*B)- JKL
            VDWHELL(3*B,3*A-2)= VDWHELL(3*B,3*A-2)- JKL
            VDWHELL(3*B-2,3*A)= VDWHELL(3*B-2,3*A)- JKL
            VDWHELL(3*A,3*B-2)= VDWHELL(3*A,3*B-2)- JKL
C
C DYADZA,DZADYA,DYBDZB,DZBDYB
C
            JKL= VIXEN2*VIXEN1*YAB*ZAB
            VDWHELL(3*A-1,3*A)= VDWHELL(3*A-1,3*A)+ JKL
            VDWHELL(3*A,3*A-1)= VDWHELL(3*A,3*A-1)+ JKL
            VDWHELL(3*B-1,3*B)= VDWHELL(3*B-1,3*B)+ JKL
            VDWHELL(3*B,3*B-1)= VDWHELL(3*B,3*B-1)+ JKL
C
C DYADZB,DZBDYA,DYBDZA,DZADYB
C
            VDWHELL(3*A-1,3*B)= VDWHELL(3*A-1,3*B)- JKL
            VDWHELL(3*B,3*A-1)= VDWHELL(3*B,3*A-1)- JKL
            VDWHELL(3*B-1,3*A)= VDWHELL(3*B-1,3*A)- JKL
            VDWHELL(3*A,3*B-1)= VDWHELL(3*A,3*B-1)- JKL
         END DO
      END DO

      CALL ANGLESECDERIVS

      RETURN

      END

   
      SUBROUTINE ANGLESECDERIVS
      USE COMMONS
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE

      DOUBLE PRECISION      VIXEN1,VIXEN2,VIXEN3,VIXEN4,VIXEN5,VIXEN6,VIXEN7,VIXEN8,VIXEN9
      DOUBLE PRECISION      U,DTHETAX1,DTHETAY1,DTHETAZ1,DTHETAX2,DTHETAY2,DTHETAZ2,DTHETAX3,DTHETAY3,DTHETAZ3,KTHETA,THETA0
      DOUBLE PRECISION      V,XP,YP,ZP,XPP,YPP,ZPP,DLAMBDAX1,DLAMBDAY1,DLAMBDAZ1,DLAMBDAX2,DLAMBDAY2,DLAMBDAZ2
      DOUBLE PRECISION      DLAMBDAX3,DLAMBDAY3,DLAMBDAZ3,DMUX2,DMUY2,DMUZ2,DMUX3,DMUY3,DMUZ3,DMUX4,DMUY4,DMUZ4
      DOUBLE PRECISION      TORSCRAP(3*NATOMS)
      DOUBLE PRECISION      MUCRAP, LAMBDACRAP, RBC, RBC2, RBC4, RBC6, XE, YE, ZE, RAB, RAB5, RAB2, RAB3
      DOUBLE PRECISION      KBOND, RBOND, RAB4, ST, CT
      DOUBLE PRECISION      DUX1,DUX2,DUX3,DUX4,DUY1,DUY2,DUY3,DUY4,DUZ1,DUZ2,DUZ3,DUZ4
      COMMON /DU/           DUX1,DUX2,DUX3,DUX4,DUY1,DUY2,DUY3,DUY4,DUZ1,DUZ2,DUZ3,DUZ4
      DOUBLE PRECISION      DVX1,DVX2,DVX3,DVX4,DVY1,DVY2,DVY3,DVY4,DVZ1,DVZ2,DVZ3,DVZ4
      COMMON /DV/           DVX1,DVX2,DVX3,DVX4,DVY1,DVY2,DVY3,DVY4,DVZ1,DVZ2,DVZ3,DVZ4
      DOUBLE PRECISION      DMX1,DMX2,DMX3,DMX4,DMY1,DMY2,DMY3,DMY4,DMZ1,DMZ2,DMZ3,DMZ4
      COMMON /DM/           DMX1,DMX2,DMX3,DMX4,DMY1,DMY2,DMY3,DMY4,DMZ1,DMZ2,DMZ3,DMZ4
      DOUBLE PRECISION      DPX1,DPX2,DPX3,DPY1,DPY2,DPY3,DPZ1,DPZ2,DPZ3
      COMMON /DP/           DPX1,DPX2,DPX3,DPY1,DPY2,DPY3,DPZ1,DPZ2,DPZ3
      DOUBLE PRECISION      DNX1,DNX2,DNX3,DNX4,DNY1,DNY2,DNY3,DNY4,DNZ1,DNZ2,DNZ3,DNZ4
      COMMON /DN/           DNX1,DNX2,DNX3,DNX4,DNY1,DNY2,DNY3,DNY4,DNZ1,DNZ2,DNZ3,DNZ4
      DOUBLE PRECISION      DQX2,DQX3,DQX4,DQY2,DQY3,DQY4,DQZ2,DQZ3,DQZ4
      COMMON /DQ/           DQX2,DQX3,DQX4,DQY2,DQY3,DQY4,DQZ2,DQZ3,DQZ4
      DOUBLE PRECISION      M,N
      COMMON /MN/           M,N

C
C ANGLE DERIVATIVES
C *****************
C
      DO I=1,ANG
         A=AA1(I)
         B=AA2(I)
         C=AA3(I) 
C
C FOR COMPARISON WITH NOTES, VIXEN1=M, VIXEN2=N, VIXEN3=DM, VIXEN4=DN
C
         KTHETA=KT(TYPE(A),TYPE(B),TYPE(C))
         THETA0=TO(TYPE(A),TYPE(B),TYPE(C))
         ST=SIN(THETA(I))
         CT=COS(THETA(I))
         U=X(B)*(X(B)-X(A)-X(C))+Y(B)*(Y(B)-Y(A)-Y(C))+Z(B)*(Z(B)-Z(A)-Z(C))+X(A)*X(C)+Y(A)*Y(C)+Z(A)*Z(C)
         VIXEN1=2.0D0*KTHETA*(THETA(I)-THETA0)/ST

         DTHETAX1=(U*(X(A)-X(B))/(R(A,B)**3*R(B,C)) - (X(C)-X(B))/(R(A,B)*R(B,C)))/ST
         DTHETAY1=(U*(Y(A)-Y(B))/(R(A,B)**3*R(B,C)) - (Y(C)-Y(B))/(R(A,B)*R(B,C)))/ST
         DTHETAZ1=(U*(Z(A)-Z(B))/(R(A,B)**3*R(B,C)) - (Z(C)-Z(B))/(R(A,B)*R(B,C)))/ST

         VIXEN5= U*(R(A,B)**2*(X(B)-X(C)) + R(B,C)**2*(X(B)-X(A))) / (R(A,B)**3*R(B,C)**3)
         VIXEN6= (2.0D0*X(B)-X(A)-X(C))/(R(A,B)*R(B,C))
         DTHETAX2=(VIXEN5-VIXEN6)/ST

         VIXEN5=U*(R(A,B)**2*(Y(B)-Y(C))+R(B,C)**2*(Y(B)-Y(A)))/(R(A,B)**3*R(B,C)**3)
         VIXEN6=(2.0D0*Y(B)-Y(A)-Y(C))/(R(A,B)*R(B,C))
         DTHETAY2=(VIXEN5-VIXEN6)/ST

         VIXEN5=U*(R(A,B)**2*(Z(B)-Z(C))+R(B,C)**2*(Z(B)-Z(A)))/(R(A,B)**3*R(B,C)**3)
         VIXEN6=(2.0D0*Z(B)-Z(A)-Z(C))/(R(A,B)*R(B,C))
         DTHETAZ2=(VIXEN5-VIXEN6)/ST

         DTHETAX3=(U*(X(C)-X(B))/(R(A,B)*R(B,C)**3)-(X(A)-X(B))/(R(A,B)*R(B,C)))/ST
         DTHETAY3=(U*(Y(C)-Y(B))/(R(A,B)*R(B,C)**3)-(Y(A)-Y(B))/(R(A,B)*R(B,C)))/ST
         DTHETAZ3=(U*(Z(C)-Z(B))/(R(A,B)*R(B,C)**3)-(Z(A)-Z(B))/(R(A,B)*R(B,C)))/ST
C
C DXADXA,DYADYA,DZADZA
C
         VIXEN2=U*(X(A)-X(B))/(R(A,B)**3*R(B,C)) - (X(C)-X(B))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAX1*(ST-(THETA(I)-THETA0)*CT)/ST**2
         VIXEN4=(U+2.0D0*(X(A)-X(B))*(X(C)-X(B))-3.0D0*U*(X(A)-X(B))**2/R(A,B)**2)/(R(A,B)**3*R(B,C))

         ANGLEHELL(3*A-2,3*A-2)= ANGLEHELL(3*A-2,3*A-2)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3 

         VIXEN2=U*(Y(A)-Y(B))/(R(A,B)**3*R(B,C)) - (Y(C)-Y(B))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAY1*(ST-(THETA(I)-THETA0)*CT)/ST**2
         VIXEN4=(U+2.0D0*(Y(A)-Y(B))*(Y(C)-Y(B))-3.0D0*U*(Y(A)-Y(B))**2/R(A,B)**2)/(R(A,B)**3*R(B,C))

         ANGLEHELL(3*A-1,3*A-1)= ANGLEHELL(3*A-1,3*A-1)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3

         VIXEN2=U*(Z(A)-Z(B))/(R(A,B)**3*R(B,C)) - (Z(C)-Z(B))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAZ1*(ST-(THETA(I)-THETA0)*CT)/ST**2
         VIXEN4=(U+2.0D0*(Z(A)-Z(B))*(Z(C)-Z(B))-3.0D0*U*(Z(A)-Z(B))**2/R(A,B)**2)/(R(A,B)**3*R(B,C))

         ANGLEHELL(3*A,3*A)= ANGLEHELL(3*A,3*A)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
C
C DXADYA,DYADXA
C
         VIXEN2=U*(X(A)-X(B))/(R(A,B)**3*R(B,C)) - (X(C)-X(B))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAY1*(ST-(THETA(I)-THETA0)*CT)/ST**2
         VIXEN4=((X(A)-X(B))*(Y(C)-Y(B))+(X(C)-X(B))*(Y(A)-Y(B))-3.0D0*U*(X(A)-X(B))*(Y(A)-Y(B))/R(A,B)**2)/(R(A,B)**3*R(B,C))
 
         ANGLEHELL(3*A-2,3*A-1)= ANGLEHELL(3*A-2,3*A-1)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
         ANGLEHELL(3*A-1,3*A-2)= ANGLEHELL(3*A-1,3*A-2)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
C
C DXADZA,DZADXA
C
         VIXEN3=2.0D0*KTHETA*DTHETAZ1*(ST-(THETA(I)-THETA0)*CT)/ST**2
         VIXEN4=((X(A)-X(B))*(Z(C)-Z(B))+(X(C)-X(B))*(Z(A)-Z(B))-3.0D0*U*(X(A)-X(B))*(Z(A)-Z(B))/R(A,B)**2)/(R(A,B)**3*R(B,C))

         ANGLEHELL(3*A-2,3*A)= ANGLEHELL(3*A-2,3*A)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
         ANGLEHELL(3*A,3*A-2)= ANGLEHELL(3*A,3*A-2)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
C
C DYADZA,DZADYA
C
         VIXEN2=U*(Y(A)-Y(B))/(R(A,B)**3*R(B,C)) - (Y(C)-Y(B))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAZ1*(ST-(THETA(I)-THETA0)*CT)/ST**2
         VIXEN4=((Y(A)-Y(B))*(Z(C)-Z(B))+(Y(C)-Y(B))*(Z(A)-Z(B))-3.0D0*U*(Y(A)-Y(B))*(Z(A)-Z(B))/R(A,B)**2)/(R(A,B)**3*R(B,C))

         ANGLEHELL(3*A-1,3*A)= ANGLEHELL(3*A-1,3*A)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
         ANGLEHELL(3*A,3*A-1)= ANGLEHELL(3*A,3*A-1)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
C
C DXADXB,DXBDXA
C
         VIXEN2=U*(X(A)-X(B))/(R(A,B)**3*R(B,C)) - (X(C)-X(B))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAX2*(ST-(THETA(I)-THETA0)*CT)/ST**2

         VIXEN5=((X(A)-X(B))*(2.0D0*X(B)-X(A)-X(C))-U)/(R(A,B)**3*R(B,C)) + 1.0D0/(R(A,B)*R(B,C))
         VIXEN6=(U*(X(A)-X(B))*(R(A,B)**3*(X(B)-X(C)) + 3.0D0*R(B,C)**2*R(A,B)*(X(B)-X(A))))/(R(A,B)**6*R(B,C)**3)
         VIXEN7=(X(C)-X(B))*(R(A,B)**2*(X(B)-X(C)) + R(B,C)**2*(X(B)-X(A)))/(R(A,B)**3*R(B,C)**3)
         VIXEN4= VIXEN5 - VIXEN6 + VIXEN7

         ANGLEHELL(3*A-2,3*B-2)= ANGLEHELL(3*A-2,3*B-2)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
         ANGLEHELL(3*B-2,3*A-2)= ANGLEHELL(3*B-2,3*A-2)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
C
C DYADYB,DYBDYA
C
         VIXEN2=U*(Y(A)-Y(B))/(R(A,B)**3*R(B,C)) - (Y(C)-Y(B))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAY2*(ST-(THETA(I)-THETA0)*CT)/ST**2

         VIXEN5=((Y(A)-Y(B))*(2.0D0*Y(B)-Y(A)-Y(C))-U)/(R(A,B)**3*R(B,C)) + 1.0D0/(R(A,B)*R(B,C))
         VIXEN6=(U*(Y(A)-Y(B))*(R(A,B)**3*(Y(B)-Y(C)) + 3.0D0*R(B,C)**2*R(A,B)*(Y(B)-Y(A))))/(R(A,B)**6*R(B,C)**3)
         VIXEN7=(Y(C)-Y(B))*(R(A,B)**2*(Y(B)-Y(C)) + R(B,C)**2*(Y(B)-Y(A)))/(R(A,B)**3*R(B,C)**3)
         VIXEN4= VIXEN5 - VIXEN6 + VIXEN7

         ANGLEHELL(3*A-1,3*B-1)= ANGLEHELL(3*A-1,3*B-1)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
         ANGLEHELL(3*B-1,3*A-1)= ANGLEHELL(3*B-1,3*A-1)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
C
C DZADZB,DZBDZA
C
         VIXEN2=U*(Z(A)-Z(B))/(R(A,B)**3*R(B,C)) - (Z(C)-Z(B))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAZ2*(ST-(THETA(I)-THETA0)*CT)/ST**2

         VIXEN5=((Z(A)-Z(B))*(2.0D0*Z(B)-Z(A)-Z(C))-U)/(R(A,B)**3*R(B,C)) + 1.0D0/(R(A,B)*R(B,C))
         VIXEN6=(U*(Z(A)-Z(B))*(R(A,B)**3*(Z(B)-Z(C)) + 3.0D0*R(B,C)**2*R(A,B)*(Z(B)-Z(A))))/(R(A,B)**6*R(B,C)**3)
         VIXEN7=(Z(C)-Z(B))*(R(A,B)**2*(Z(B)-Z(C)) + R(B,C)**2*(Z(B)-Z(A)))/(R(A,B)**3*R(B,C)**3)
         VIXEN4= VIXEN5 - VIXEN6 + VIXEN7

         ANGLEHELL(3*A,3*B)= ANGLEHELL(3*A,3*B)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
         ANGLEHELL(3*B,3*A)= ANGLEHELL(3*B,3*A)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
C
C DXADYB,DYBDXA
C
         VIXEN2=U*(X(A)-X(B))/(R(A,B)**3*R(B,C)) - (X(C)-X(B))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAY2*(ST-(THETA(I)-THETA0)*CT)/ST**2

         VIXEN5=(X(A)-X(B))*(2.0D0*Y(B)-Y(A)-Y(C))/(R(A,B)**3*R(B,C)) 
         VIXEN6=U*(X(A)-X(B))*(R(A,B)**3*(Y(B)-Y(C)) + 3.0D0*R(A,B)*R(B,C)**2*(Y(B)-Y(A)))/(R(A,B)**6*R(B,C)**3)
         VIXEN4= VIXEN5 - VIXEN6 + (X(C)-X(B))*((Y(B)-Y(C))/(R(A,B)*R(B,C)**3) + (Y(B)-Y(A))/(R(A,B)**3*R(B,C)))

         ANGLEHELL(3*A-2,3*B-1)= ANGLEHELL(3*A-2,3*B-1)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
         ANGLEHELL(3*B-1,3*A-2)= ANGLEHELL(3*B-1,3*A-2)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
C
C DXADZB,DZBDXA
C
         VIXEN2=U*(X(A)-X(B))/(R(A,B)**3*R(B,C)) - (X(C)-X(B))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAZ2*(ST-(THETA(I)-THETA0)*CT)/ST**2

         VIXEN5=(X(A)-X(B))*(2.0D0*Z(B)-Z(A)-Z(C))/(R(A,B)**3*R(B,C)) 
         VIXEN6=U*(X(A)-X(B))*(R(A,B)**3*(Z(B)-Z(C))+3.0D0*R(A,B)*R(B,C)**2*(Z(B)-Z(A)))/(R(A,B)**6*R(B,C)**3)
         VIXEN4= VIXEN5 - VIXEN6 + (X(C)-X(B))*((Z(B)-Z(C))/(R(A,B)*R(B,C)**3) + (Z(B)-Z(A))/(R(A,B)**3*R(B,C)))

         ANGLEHELL(3*A-2,3*B)= ANGLEHELL(3*A-2,3*B)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
         ANGLEHELL(3*B,3*A-2)= ANGLEHELL(3*B,3*A-2)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
C
C DYADXB,DXBDYA
C
         VIXEN2=U*(Y(A)-Y(B))/(R(A,B)**3*R(B,C)) - (Y(C)-Y(B))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAX2*(ST-(THETA(I)-THETA0)*CT)/ST**2

         VIXEN5=(Y(A)-Y(B))*(2.0D0*X(B)-X(A)-X(C))/(R(A,B)**3*R(B,C))
         VIXEN6=U*(Y(A)-Y(B))*(R(A,B)**3*(X(B)-X(C))+3.0D0*R(A,B)*R(B,C)**2*(X(B)-X(A)))/(R(A,B)**6*R(B,C)**3)
         VIXEN4= VIXEN5 - VIXEN6 + (Y(C)-Y(B))*((X(B)-X(C))/(R(A,B)*R(B,C)**3) + (X(B)-X(A))/(R(A,B)**3*R(B,C)))

         ANGLEHELL(3*A-1,3*B-2)= ANGLEHELL(3*A-1,3*B-2)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
         ANGLEHELL(3*B-2,3*A-1)= ANGLEHELL(3*B-2,3*A-1)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
C
C DYADZB,DZBDYA
C
         VIXEN2=U*(Y(A)-Y(B))/(R(A,B)**3*R(B,C)) - (Y(C)-Y(B))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAZ2*(ST-(THETA(I)-THETA0)*CT)/ST**2

         VIXEN5=(Y(A)-Y(B))*(2.0D0*Z(B)-Z(A)-Z(C))/(R(A,B)**3*R(B,C))
         VIXEN6=U*(Y(A)-Y(B))*(R(A,B)**3*(Z(B)-Z(C))+3.0D0*R(A,B)*R(B,C)**2*(Z(B)-Z(A)))/(R(A,B)**6*R(B,C)**3)
         VIXEN4= VIXEN5 - VIXEN6 + (Y(C)-Y(B))*((Z(B)-Z(C))/(R(A,B)*R(B,C)**3) + (Z(B)-Z(A))/(R(A,B)**3*R(B,C)))

         ANGLEHELL(3*A-1,3*B)= ANGLEHELL(3*A-1,3*B)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
         ANGLEHELL(3*B,3*A-1)= ANGLEHELL(3*B,3*A-1)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
C
C DZADXB,DXBDZA
C
         VIXEN2=U*(Z(A)-Z(B))/(R(A,B)**3*R(B,C)) - (Z(C)-Z(B))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAX2*(ST-(THETA(I)-THETA0)*CT)/ST**2

         VIXEN5=(Z(A)-Z(B))*(2.0D0*X(B)-X(A)-X(C))/(R(A,B)**3*R(B,C))
         VIXEN6=U*(Z(A)-Z(B))*(R(A,B)**3*(X(B)-X(C))+3.0D0*R(A,B)*R(B,C)**2*(X(B)-X(A)))/(R(A,B)**6*R(B,C)**3)
         VIXEN4= VIXEN5 - VIXEN6 + (Z(C)-Z(B))*((X(B)-X(C))/(R(A,B)*R(B,C)**3) + (X(B)-X(A))/(R(A,B)**3*R(B,C)))

         ANGLEHELL(3*A,3*B-2)= ANGLEHELL(3*A,3*B-2)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
         ANGLEHELL(3*B-2,3*A)= ANGLEHELL(3*B-2,3*A)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
C
C DZADYB,DYBDZA
C
         VIXEN2=U*(Z(A)-Z(B))/(R(A,B)**3*R(B,C)) - (Z(C)-Z(B))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAY2*(ST-(THETA(I)-THETA0)*CT)/ST**2

         VIXEN5=(Z(A)-Z(B))*(2.0D0*Y(B)-Y(A)-Y(C))/(R(A,B)**3*R(B,C))
         VIXEN6=U*(Z(A)-Z(B))*(R(A,B)**3*(Y(B)-Y(C))+3.0D0*R(A,B)*R(B,C)**2*(Y(B)-Y(A)))/(R(A,B)**6*R(B,C)**3)
         VIXEN4= VIXEN5 - VIXEN6 + (Z(C)-Z(B))*((Y(B)-Y(C))/(R(A,B)*R(B,C)**3) + (Y(B)-Y(A))/(R(A,B)**3*R(B,C)))

         ANGLEHELL(3*A,3*B-1)= ANGLEHELL(3*A,3*B-1)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
         ANGLEHELL(3*B-1,3*A)= ANGLEHELL(3*B-1,3*A)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
C
C DXADXC,DXCDXA
C
         VIXEN2=U*(X(A)-X(B))/(R(A,B)**3*R(B,C)) - (X(C)-X(B))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAX3*(ST-(THETA(I)-THETA0)*CT)/ST**2

         VIXEN5=(X(A)-X(B))**2/(R(A,B)**3*R(B,C)) + (X(C)-X(B))**2/(R(A,B)*R(B,C)**3) - 1.0D0/(R(A,B)*R(B,C))
         VIXEN4= VIXEN5 - U*(X(A)-X(B))*(X(C)-X(B))/(R(A,B)**3*R(B,C)**3)

         ANGLEHELL(3*A-2,3*C-2)= ANGLEHELL(3*A-2,3*C-2)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
         ANGLEHELL(3*C-2,3*A-2)= ANGLEHELL(3*C-2,3*A-2)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
C
C DYADYC,DYCDYA
C
         VIXEN2=U*(Y(A)-Y(B))/(R(A,B)**3*R(B,C)) - (Y(C)-Y(B))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAY3*(ST-(THETA(I)-THETA0)*CT)/ST**2

         VIXEN5=(Y(A)-Y(B))**2/(R(A,B)**3*R(B,C)) + (Y(C)-Y(B))**2/(R(A,B)*R(B,C)**3) - 1.0D0/(R(A,B)*R(B,C))
         VIXEN4= VIXEN5 - U*(Y(A)-Y(B))*(Y(C)-Y(B))/(R(A,B)**3*R(B,C)**3)

         ANGLEHELL(3*A-1,3*C-1)= ANGLEHELL(3*A-1,3*C-1)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
         ANGLEHELL(3*C-1,3*A-1)= ANGLEHELL(3*C-1,3*A-1)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
C
C DZADZC,DZCDZA
C
         VIXEN2=U*(Z(A)-Z(B))/(R(A,B)**3*R(B,C)) - (Z(C)-Z(B))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAZ3*(ST-(THETA(I)-THETA0)*CT)/ST**2

         VIXEN5=(Z(A)-Z(B))**2/(R(A,B)**3*R(B,C)) + (Z(C)-Z(B))**2/(R(A,B)*R(B,C)**3) - 1.0D0/(R(A,B)*R(B,C))
         VIXEN4= VIXEN5 - U*(Z(A)-Z(B))*(Z(C)-Z(B))/(R(A,B)**3*R(B,C)**3)

         ANGLEHELL(3*A,3*C)= ANGLEHELL(3*A,3*C)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
         ANGLEHELL(3*C,3*A)= ANGLEHELL(3*C,3*A)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
C
C DXADYC,DYCDXA
C
         VIXEN2=U*(X(A)-X(B))/(R(A,B)**3*R(B,C)) - (X(C)-X(B))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAY3*(ST-(THETA(I)-THETA0)*CT)/ST**2

         VIXEN5=(X(A)-X(B))*(Y(A)-Y(B))/(R(A,B)**3*R(B,C)) - U*(X(A)-X(B))*(Y(C)-Y(B))/(R(A,B)**3*R(B,C)**3) 
         VIXEN4= VIXEN5 + (X(C)-X(B))*(Y(C)-Y(B))/(R(A,B)*R(B,C)**3)

         ANGLEHELL(3*A-2,3*C-1)= ANGLEHELL(3*A-2,3*C-1)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
         ANGLEHELL(3*C-1,3*A-2)= ANGLEHELL(3*C-1,3*A-2)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
C
C DXADZC,DZCDXA
C
         VIXEN2=U*(X(A)-X(B))/(R(A,B)**3*R(B,C)) - (X(C)-X(B))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAZ3*(ST-(THETA(I)-THETA0)*CT)/ST**2

         VIXEN5=(X(A)-X(B))*(Z(A)-Z(B))/(R(A,B)**3*R(B,C)) - U*(X(A)-X(B))*(Z(C)-Z(B))/(R(A,B)**3*R(B,C)**3) 
         VIXEN4= VIXEN5 + (X(C)-X(B))*(Z(C)-Z(B))/(R(A,B)*R(B,C)**3)

         ANGLEHELL(3*A-2,3*C)= ANGLEHELL(3*A-2,3*C)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
         ANGLEHELL(3*C,3*A-2)= ANGLEHELL(3*C,3*A-2)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
C
C DYADZC,DZCDYA
C
         VIXEN2=U*(Y(A)-Y(B))/(R(A,B)**3*R(B,C)) - (Y(C)-Y(B))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAZ3*(ST-(THETA(I)-THETA0)*CT)/ST**2

         VIXEN5=(Y(A)-Y(B))*(Z(A)-Z(B))/(R(A,B)**3*R(B,C)) - U*(Y(A)-Y(B))*(Z(C)-Z(B))/(R(A,B)**3*R(B,C)**3) 
         VIXEN4= VIXEN5 + (Y(C)-Y(B))*(Z(C)-Z(B))/(R(A,B)*R(B,C)**3)

         ANGLEHELL(3*A-1,3*C)= ANGLEHELL(3*A-1,3*C)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
         ANGLEHELL(3*C,3*A-1)= ANGLEHELL(3*C,3*A-1)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
C
C DYADXC,DXCDYA
C
         VIXEN2=U*(Y(A)-Y(B))/(R(A,B)**3*R(B,C)) - (Y(C)-Y(B))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAX3*(ST-(THETA(I)-THETA0)*CT)/ST**2

         VIXEN5=(Y(A)-Y(B))*(X(A)-X(B))/(R(A,B)**3*R(B,C)) - U*(Y(A)-Y(B))*(X(C)-X(B))/(R(A,B)**3*R(B,C)**3) 
         VIXEN4= VIXEN5 + (Y(C)-Y(B))*(X(C)-X(B))/(R(A,B)*R(B,C)**3)

         ANGLEHELL(3*A-1,3*C-2)= ANGLEHELL(3*A-1,3*C-2)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
         ANGLEHELL(3*C-2,3*A-1)= ANGLEHELL(3*C-2,3*A-1)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
C
C DZADXC,DXCDZA
C
         VIXEN2=U*(Z(A)-Z(B))/(R(A,B)**3*R(B,C)) - (Z(C)-Z(B))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAX3*(ST-(THETA(I)-THETA0)*CT)/ST**2

         VIXEN5=(Z(A)-Z(B))*(X(A)-X(B))/(R(A,B)**3*R(B,C)) - U*(Z(A)-Z(B))*(X(C)-X(B))/(R(A,B)**3*R(B,C)**3) 
         VIXEN4= VIXEN5 + (Z(C)-Z(B))*(X(C)-X(B))/(R(A,B)*R(B,C)**3)

         ANGLEHELL(3*A,3*C-2)= ANGLEHELL(3*A,3*C-2)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
         ANGLEHELL(3*C-2,3*A)= ANGLEHELL(3*C-2,3*A)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
C
C DZADYC,DYCDZA
C
         VIXEN2=U*(Z(A)-Z(B))/(R(A,B)**3*R(B,C)) - (Z(C)-Z(B))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAY3*(ST-(THETA(I)-THETA0)*CT)/ST**2

         VIXEN5=(Z(A)-Z(B))*(Y(A)-Y(B))/(R(A,B)**3*R(B,C)) - U*(Z(A)-Z(B))*(Y(C)-Y(B))/(R(A,B)**3*R(B,C)**3) 
         VIXEN4= VIXEN5 + (Z(C)-Z(B))*(Y(C)-Y(B))/(R(A,B)*R(B,C)**3)

         ANGLEHELL(3*A,3*C-1)= ANGLEHELL(3*A,3*C-1)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
         ANGLEHELL(3*C-1,3*A)= ANGLEHELL(3*C-1,3*A)+ VIXEN1*VIXEN4+VIXEN2*VIXEN3
C
C DXCDXC
C
         VIXEN2=U*(X(C)-X(B))/(R(A,B)*R(B,C)**3) - (X(A)-X(B))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAX3*(ST-(THETA(I)-THETA0)*CT)/ST**2
         VIXEN4=(2.0D0*(X(A)-X(B))*(X(C)-X(B)) + U - 3.0D0*U*(X(C)-X(B))**2/R(B,C)**2)/(R(A,B)*R(B,C)**3)

         ANGLEHELL(3*C-2,3*C-2)= ANGLEHELL(3*C-2,3*C-2) + VIXEN1*VIXEN4+VIXEN2*VIXEN3
C
C DYCDYC
C
         VIXEN2=U*(Y(C)-Y(B))/(R(A,B)*R(B,C)**3) - (Y(A)-Y(B))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAY3*(ST-(THETA(I)-THETA0)*CT)/ST**2
         VIXEN4=(2.0D0*(Y(A)-Y(B))*(Y(C)-Y(B)) + U - 3.0D0*U*(Y(C)-Y(B))**2/R(B,C)**2)/(R(A,B)*R(B,C)**3)

         ANGLEHELL(3*C-1,3*C-1)= ANGLEHELL(3*C-1,3*C-1) + VIXEN1*VIXEN4+VIXEN2*VIXEN3
C
C DZCDZC
C
         VIXEN2=U*(Z(C)-Z(B))/(R(A,B)*R(B,C)**3) - (Z(A)-Z(B))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAZ3*(ST-(THETA(I)-THETA0)*CT)/ST**2
         VIXEN4=(2.0D0*(Z(A)-Z(B))*(Z(C)-Z(B)) + U - 3.0D0*U*(Z(C)-Z(B))**2/R(B,C)**2)/(R(A,B)*R(B,C)**3)

         ANGLEHELL(3*C,3*C)= ANGLEHELL(3*C,3*C) + VIXEN1*VIXEN4+VIXEN2*VIXEN3
C
C DXCDYC,DYCDXC
C
         VIXEN2=U*(X(C)-X(B))/(R(A,B)*R(B,C)**3) - (X(A)-X(B))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAY3*(ST-(THETA(I)-THETA0)*CT)/ST**2
         VIXEN4=((X(C)-X(B))*(Y(A)-Y(B))+(X(A)-X(B))*(Y(C)-Y(B)) - 3.0D0*U*(X(C)-X(B))*(Y(C)-Y(B))/R(B,C)**2)/(R(A,B)*R(B,C)**3)

         ANGLEHELL(3*C-2,3*C-1)=ANGLEHELL(3*C-2,3*C-1) + VIXEN1*VIXEN4+VIXEN2*VIXEN3
         ANGLEHELL(3*C-1,3*C-2)=ANGLEHELL(3*C-1,3*C-2) + VIXEN1*VIXEN4+VIXEN2*VIXEN3 
C
C DXCDZC,DZCDXC
C
         VIXEN2=U*(X(C)-X(B))/(R(A,B)*R(B,C)**3) - (X(A)-X(B))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAZ3*(ST-(THETA(I)-THETA0)*CT)/ST**2
         VIXEN4=((X(C)-X(B))*(Z(A)-Z(B))+(X(A)-X(B))*(Z(C)-Z(B)) - 3.0D0*U*(X(C)-X(B))*(Z(C)-Z(B))/R(B,C)**2)/(R(A,B)*R(B,C)**3)

         ANGLEHELL(3*C-2,3*C)=ANGLEHELL(3*C-2,3*C) + VIXEN1*VIXEN4+VIXEN2*VIXEN3
         ANGLEHELL(3*C,3*C-2)=ANGLEHELL(3*C,3*C-2) + VIXEN1*VIXEN4+VIXEN2*VIXEN3
C
C DYCDZC,DZCDYC
C
         VIXEN2=U*(Y(C)-Y(B))/(R(A,B)*R(B,C)**3) - (Y(A)-Y(B))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAZ3*(ST-(THETA(I)-THETA0)*CT)/ST**2
         VIXEN4=((Y(C)-Y(B))*(Z(A)-Z(B))+(Y(A)-Y(B))*(Z(C)-Z(B)) - 3.0D0*U*(Y(C)-Y(B))*(Z(C)-Z(B))/R(B,C)**2)/(R(A,B)*R(B,C)**3)

         ANGLEHELL(3*C-1,3*C)=ANGLEHELL(3*C-1,3*C) + VIXEN1*VIXEN4+VIXEN2*VIXEN3
         ANGLEHELL(3*C,3*C-1)=ANGLEHELL(3*C,3*C-1) + VIXEN1*VIXEN4+VIXEN2*VIXEN3
C
C DXCDXB,DXBDXC
C
         VIXEN2=U*(X(C)-X(B))/(R(A,B)*R(B,C)**3) - (X(A)-X(B))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAX2*(ST-(THETA(I)-THETA0)*CT)/ST**2

         VIXEN5=((X(C)-X(B))*(2.0D0*X(B)-X(C)-X(A))-U)/(R(A,B)*R(B,C)**3) + 1.0D0/(R(A,B)*R(B,C))
         VIXEN6=U*(X(C)-X(B))*(3.0D0*R(A,B)**2*(X(B)-X(C)) + R(B,C)**2*(X(B)-X(A)))/(R(A,B)**3*R(B,C)**5)
         VIXEN4= VIXEN5 - VIXEN6 + (X(A)-X(B))*(R(A,B)**2*(X(B)-X(C)) + R(B,C)**2*(X(B)-X(A)))/(R(A,B)**3*R(B,C)**3)

         ANGLEHELL(3*C-2,3*B-2)= ANGLEHELL(3*C-2,3*B-2) + VIXEN1*VIXEN4+VIXEN2*VIXEN3
         ANGLEHELL(3*B-2,3*C-2)= ANGLEHELL(3*B-2,3*C-2) + VIXEN1*VIXEN4+VIXEN2*VIXEN3
C
C DYCDYB,DYBDYC
C
         VIXEN2=U*(Y(C)-Y(B))/(R(A,B)*R(B,C)**3) - (Y(A)-Y(B))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAY2*(ST-(THETA(I)-THETA0)*CT)/ST**2

         VIXEN5=((Y(C)-Y(B))*(2.0D0*Y(B)-Y(C)-Y(A))-U)/(R(A,B)*R(B,C)**3) + 1.0D0/(R(A,B)*R(B,C))
         VIXEN6=U*(Y(C)-Y(B))*(3.0D0*R(A,B)**2*(Y(B)-Y(C)) + R(B,C)**2*(Y(B)-Y(A)))/(R(A,B)**3*R(B,C)**5)
         VIXEN4= VIXEN5 - VIXEN6 + (Y(A)-Y(B))*(R(A,B)**2*(Y(B)-Y(C)) + R(B,C)**2*(Y(B)-Y(A)))/(R(A,B)**3*R(B,C)**3)

         ANGLEHELL(3*C-1,3*B-1)= ANGLEHELL(3*C-1,3*B-1) + VIXEN1*VIXEN4+VIXEN2*VIXEN3
         ANGLEHELL(3*B-1,3*C-1)= ANGLEHELL(3*B-1,3*C-1) + VIXEN1*VIXEN4+VIXEN2*VIXEN3
C
C DZCDZB,DZBDZC
C
         VIXEN2=U*(Z(C)-Z(B))/(R(A,B)*R(B,C)**3) - (Z(A)-Z(B))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAZ2*(ST-(THETA(I)-THETA0)*CT)/ST**2

         VIXEN5=((Z(C)-Z(B))*(2.0D0*Z(B)-Z(C)-Z(A))-U)/(R(A,B)*R(B,C)**3) + 1.0D0/(R(A,B)*R(B,C))
         VIXEN6=U*(Z(C)-Z(B))*(3.0D0*R(A,B)**2*(Z(B)-Z(C)) + R(B,C)**2*(Z(B)-Z(A)))/(R(A,B)**3*R(B,C)**5)
         VIXEN4= VIXEN5 - VIXEN6 + (Z(A)-Z(B))*(R(A,B)**2*(Z(B)-Z(C)) + R(B,C)**2*(Z(B)-Z(A)))/(R(A,B)**3*R(B,C)**3)

         ANGLEHELL(3*C,3*B)= ANGLEHELL(3*C,3*B) + VIXEN1*VIXEN4+VIXEN2*VIXEN3
         ANGLEHELL(3*B,3*C)= ANGLEHELL(3*B,3*C) + VIXEN1*VIXEN4+VIXEN2*VIXEN3
C
C DXCDYB,DYBDXC
C
         VIXEN2=U*(X(C)-X(B))/(R(A,B)*R(B,C)**3) - (X(A)-X(B))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAY2*(ST-(THETA(I)-THETA0)*CT)/ST**2 

         VIXEN5=(X(C)-X(B))*(2.0D0*Y(B)-Y(C)-Y(A))/(R(A,B)*R(B,C)**3) 
         VIXEN6=U*(X(C)-X(B))*(3.0D0*R(A,B)**2*(Y(B)-Y(C)) + R(B,C)**2*(Y(B)-Y(A)))/(R(A,B)**3*R(B,C)**5)
         VIXEN4= VIXEN5 - VIXEN6 + (X(A)-X(B))*((Y(B)-Y(C))/R(B,C)**2 + (Y(B)-Y(A))/R(A,B)**2)/(R(A,B)*R(B,C))

         ANGLEHELL(3*C-2,3*B-1)= ANGLEHELL(3*C-2,3*B-1) + VIXEN1*VIXEN4+VIXEN2*VIXEN3
         ANGLEHELL(3*B-1,3*C-2)= ANGLEHELL(3*B-1,3*C-2) + VIXEN1*VIXEN4+VIXEN2*VIXEN3
C
C DXCDZB,DZBDXC
C
         VIXEN2=U*(X(C)-X(B))/(R(A,B)*R(B,C)**3) - (X(A)-X(B))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAZ2*(ST-(THETA(I)-THETA0)*CT)/ST**2

         VIXEN5=(X(C)-X(B))*(2.0D0*Z(B)-Z(C)-Z(A))/(R(A,B)*R(B,C)**3) 
         VIXEN6=U*(X(C)-X(B))*(3.0D0*R(A,B)**2*(Z(B)-Z(C)) + R(B,C)**2*(Z(B)-Z(A)))/(R(A,B)**3*R(B,C)**5)
         VIXEN4= VIXEN5 - VIXEN6 + (X(A)-X(B))*((Z(B)-Z(C))/R(B,C)**2 + (Z(B)-Z(A))/R(A,B)**2)/(R(A,B)*R(B,C))

         ANGLEHELL(3*C-2,3*B)= ANGLEHELL(3*C-2,3*B) + VIXEN1*VIXEN4+VIXEN2*VIXEN3
         ANGLEHELL(3*B,3*C-2)= ANGLEHELL(3*B,3*C-2) + VIXEN1*VIXEN4+VIXEN2*VIXEN3
C
C DYCDZB,DZBDYC
C
         VIXEN2=U*(Y(C)-Y(B))/(R(A,B)*R(B,C)**3) - (Y(A)-Y(B))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAZ2*(ST-(THETA(I)-THETA0)*CT)/ST**2

         VIXEN5=(Y(C)-Y(B))*(2.0D0*Z(B)-Z(C)-Z(A))/(R(A,B)*R(B,C)**3) 
         VIXEN6=U*(Y(C)-Y(B))*(3.0D0*R(A,B)**2*(Z(B)-Z(C)) + R(B,C)**2*(Z(B)-Z(A)))/(R(A,B)**3*R(B,C)**5)
         VIXEN4= VIXEN5 - VIXEN6 + (Y(A)-Y(B))*((Z(B)-Z(C))/R(B,C)**2 + (Z(B)-Z(A))/R(A,B)**2)/(R(A,B)*R(B,C))

         ANGLEHELL(3*C-1,3*B)= ANGLEHELL(3*C-1,3*B) + VIXEN1*VIXEN4+VIXEN2*VIXEN3
         ANGLEHELL(3*B,3*C-1)= ANGLEHELL(3*B,3*C-1) + VIXEN1*VIXEN4+VIXEN2*VIXEN3
C
C DYCDXB,DXBDYC
C
         VIXEN2=U*(Y(C)-Y(B))/(R(A,B)*R(B,C)**3) - (Y(A)-Y(B))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAX2*(ST-(THETA(I)-THETA0)*CT)/ST**2

         VIXEN5=(Y(C)-Y(B))*(2.0D0*X(B)-X(C)-X(A))/(R(A,B)*R(B,C)**3)
         VIXEN6=U*(Y(C)-Y(B))*(3.0D0*R(A,B)**2*(X(B)-X(C)) + R(B,C)**2*(X(B)-X(A)))/(R(A,B)**3*R(B,C)**5)
         VIXEN4= VIXEN5 - VIXEN6 + (Y(A)-Y(B))*((X(B)-X(C))/R(B,C)**2 + (X(B)-X(A))/R(A,B)**2)/(R(A,B)*R(B,C))

         ANGLEHELL(3*C-1,3*B-2)= ANGLEHELL(3*C-1,3*B-2) + VIXEN1*VIXEN4+VIXEN2*VIXEN3
         ANGLEHELL(3*B-2,3*C-1)= ANGLEHELL(3*B-2,3*C-1) + VIXEN1*VIXEN4+VIXEN2*VIXEN3 
C
C DZCDXB,DXBDZC 
C
         VIXEN2=U*(Z(C)-Z(B))/(R(A,B)*R(B,C)**3) - (Z(A)-Z(B))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAX2*(ST-(THETA(I)-THETA0)*CT)/ST**2

         VIXEN5=(Z(C)-Z(B))*(2.0D0*X(B)-X(C)-X(A))/(R(A,B)*R(B,C)**3)
         VIXEN6=U*(Z(C)-Z(B))*(3.0D0*R(A,B)**2*(X(B)-X(C)) + R(B,C)**2*(X(B)-X(A)))/(R(A,B)**3*R(B,C)**5)
         VIXEN4= VIXEN5 - VIXEN6 + (Z(A)-Z(B))*((X(B)-X(C))/R(B,C)**2 + (X(B)-X(A))/R(A,B)**2)/(R(A,B)*R(B,C))

         ANGLEHELL(3*C,3*B-2)= ANGLEHELL(3*C,3*B-2) + VIXEN1*VIXEN4+VIXEN2*VIXEN3
         ANGLEHELL(3*B-2,3*C)= ANGLEHELL(3*B-2,3*C) + VIXEN1*VIXEN4+VIXEN2*VIXEN3
C
C DZCDYB,DYBDZC 
C
         VIXEN2=U*(Z(C)-Z(B))/(R(A,B)*R(B,C)**3) - (Z(A)-Z(B))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAY2*(ST-(THETA(I)-THETA0)*CT)/ST**2

         VIXEN5=(Z(C)-Z(B))*(2.0D0*Y(B)-Y(C)-Y(A))/(R(A,B)*R(B,C)**3)
         VIXEN6=U*(Z(C)-Z(B))*(3.0D0*R(A,B)**2*(Y(B)-Y(C)) + R(B,C)**2*(Y(B)-Y(A)))/(R(A,B)**3*R(B,C)**5)
         VIXEN4= VIXEN5 - VIXEN6 + (Z(A)-Z(B))*((Y(B)-Y(C))/R(B,C)**2 + (Y(B)-Y(A))/R(A,B)**2)/(R(A,B)*R(B,C))

         ANGLEHELL(3*C,3*B-1)= ANGLEHELL(3*C,3*B-1) + VIXEN1*VIXEN4+VIXEN2*VIXEN3
         ANGLEHELL(3*B-1,3*C)= ANGLEHELL(3*B-1,3*C) + VIXEN1*VIXEN4+VIXEN2*VIXEN3
C
C DXBDXB
C
         VIXEN2=U*(R(A,B)**2*(X(B)-X(C)) + R(B,C)**2*(X(B)-X(A)))/(R(A,B)**3*R(B,C)**3) - (2.0D0*X(B)-X(C)-X(A))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAX2*(ST-(THETA(I)-THETA0)*CT)/ST**2

         VIXEN5=U*(4.0D0*(X(B)-X(C))*(X(B)-X(A))+R(A,B)**2+R(B,C)**2) + 
     1          (2.0D0*X(B)-X(C)-X(A))*((X(B)-X(C))*R(A,B)**2+(X(B)-X(A))*R(B,C)**2)
         VIXEN6=3.0D0*R(A,B)*R(B,C)*(R(A,B)**2*(X(B)-X(C)) + R(B,C)**2*(X(B)-X(A)))
         VIXEN7=U*(R(A,B)**2*(X(B)-X(C)) + R(B,C)**2*(X(B)-X(A)))
         VIXEN8= (VIXEN5*R(A,B)**3*R(B,C)**3 - VIXEN7*VIXEN6)/(R(A,B)**6*R(B,C)**6) 
         VIXEN9= (2.0D0*X(B)-X(C)-X(A))*(R(A,B)**2*(X(B)-X(C)) + R(B,C)**2*(X(B)-X(A)))/(R(A,B)**3*R(B,C)**3)
         VIXEN4= VIXEN8 - 2.0D0/(R(A,B)*R(B,C)) + VIXEN9

         ANGLEHELL(3*B-2,3*B-2)= ANGLEHELL(3*B-2,3*B-2) + VIXEN1*VIXEN4+VIXEN2*VIXEN3
C
C DYBDYB
C
         VIXEN2=U*(R(A,B)**2*(Y(B)-Y(C)) + R(B,C)**2*(Y(B)-Y(A)))/(R(A,B)**3*R(B,C)**3) - (2.0D0*Y(B)-Y(C)-Y(A))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAY2*(ST-(THETA(I)-THETA0)*CT)/ST**2

         VIXEN5=U*(4.0D0*(Y(B)-Y(C))*(Y(B)-Y(A))+R(A,B)**2+R(B,C)**2) + 
     1          (2.0D0*Y(B)-Y(C)-Y(A))*((Y(B)-Y(C))*R(A,B)**2+(Y(B)-Y(A))*R(B,C)**2)
         VIXEN6=3.0D0*R(A,B)*R(B,C)*(R(A,B)**2*(Y(B)-Y(C)) + R(B,C)**2*(Y(B)-Y(A)))
         VIXEN7=U*(R(A,B)**2*(Y(B)-Y(C)) + R(B,C)**2*(Y(B)-Y(A)))
         VIXEN8= (VIXEN5*R(A,B)**3*R(B,C)**3 - VIXEN7*VIXEN6)/(R(A,B)**6*R(B,C)**6)
         VIXEN9= (2.0D0*Y(B)-Y(C)-Y(A))*(R(A,B)**2*(Y(B)-Y(C)) + R(B,C)**2*(Y(B)-Y(A)))/(R(A,B)**3*R(B,C)**3)
         VIXEN4= VIXEN8 - 2.0D0/(R(A,B)*R(B,C)) + VIXEN9

         ANGLEHELL(3*B-1,3*B-1)= ANGLEHELL(3*B-1,3*B-1) + VIXEN1*VIXEN4+VIXEN2*VIXEN3
C
C DZBDZB
C
         VIXEN2=U*(R(A,B)**2*(Z(B)-Z(C)) + R(B,C)**2*(Z(B)-Z(A)))/(R(A,B)**3*R(B,C)**3) - (2.0D0*Z(B)-Z(C)-Z(A))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAZ2*(ST-(THETA(I)-THETA0)*CT)/ST**2

         VIXEN5=U*(4.0D0*(Z(B)-Z(C))*(Z(B)-Z(A))+R(A,B)**2+R(B,C)**2) + 
     1          (2.0D0*Z(B)-Z(C)-Z(A))*((Z(B)-Z(C))*R(A,B)**2+(Z(B)-Z(A))*R(B,C)**2)
         VIXEN6=3.0D0*R(A,B)*R(B,C)*(R(A,B)**2*(Z(B)-Z(C)) + R(B,C)**2*(Z(B)-Z(A)))
         VIXEN7=U*(R(A,B)**2*(Z(B)-Z(C)) + R(B,C)**2*(Z(B)-Z(A)))
         VIXEN8= (VIXEN5*R(A,B)**3*R(B,C)**3 - VIXEN7*VIXEN6)/(R(A,B)**6*R(B,C)**6)
         VIXEN9= (2.0D0*Z(B)-Z(C)-Z(A))*(R(A,B)**2*(Z(B)-Z(C)) + R(B,C)**2*(Z(B)-Z(A)))/(R(A,B)**3*R(B,C)**3)
         VIXEN4= VIXEN8 - 2.0D0/(R(A,B)*R(B,C)) + VIXEN9

         ANGLEHELL(3*B,3*B)= ANGLEHELL(3*B,3*B) + VIXEN1*VIXEN4+VIXEN2*VIXEN3
C
C DXBDYB,DYBDXB
C
         VIXEN2=U*(R(A,B)**2*(X(B)-X(C)) + R(B,C)**2*(X(B)-X(A)))/(R(A,B)**3*R(B,C)**3) - (2.0D0*X(B)-X(C)-X(A))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAY2*(ST-(THETA(I)-THETA0)*CT)/ST**2

         VIXEN5=U*(2.0D0*(Y(B)-Y(A))*(X(B)-X(C))+2.0D0*(Y(B)-Y(C))*(X(B)-X(A)))
         VIXEN5= VIXEN5 + (2.0D0*Y(B)-Y(C)-Y(A))*(R(A,B)**2*(X(B)-X(C))+R(B,C)**2*(X(B)-X(A)))
C
C ABOVE IS JUST A TRICK AS OTHERWISE GET 'UNBALANCED PARENTHESES' ERROR MESSAGE...
C
         VIXEN6=3.0D0*R(A,B)*R(B,C)*(R(A,B)**2*(Y(B)-Y(C)) + R(B,C)**2*(Y(B)-Y(A)))
         VIXEN7=U*(R(A,B)**2*(X(B)-X(C)) + R(B,C)**2*(X(B)-X(A)))
         VIXEN8= (VIXEN5*R(A,B)**3*R(B,C)**3 - VIXEN6*VIXEN7)/(R(A,B)**6*R(B,C)**6)
         VIXEN9=(2.0D0*X(B)-X(C)-X(A))*((Y(B)-Y(C))/R(B,C)**2 + (Y(B)-Y(A))/R(A,B)**2)/(R(A,B)*R(B,C))
         VIXEN4= VIXEN8 + VIXEN9

         ANGLEHELL(3*B-2,3*B-1)= ANGLEHELL(3*B-2,3*B-1) + VIXEN1*VIXEN4+VIXEN2*VIXEN3
         ANGLEHELL(3*B-1,3*B-2)= ANGLEHELL(3*B-1,3*B-2) + VIXEN1*VIXEN4+VIXEN2*VIXEN3
C
C DXBDZB,DZBDXB
C
         VIXEN2=U*(R(A,B)**2*(X(B)-X(C)) + R(B,C)**2*(X(B)-X(A)))/(R(A,B)**3*R(B,C)**3) - (2.0D0*X(B)-X(C)-X(A))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAZ2*(ST-(THETA(I)-THETA0)*CT)/ST**2

         VIXEN5=U*(2.0D0*(Z(B)-Z(A))*(X(B)-X(C))+2.0D0*(Z(B)-Z(C))*(X(B)-X(A)))
         VIXEN5= VIXEN5 + (2.0D0*Z(B)-Z(C)-Z(A))*(R(A,B)**2*(X(B)-X(C))+R(B,C)**2*(X(B)-X(A)))
C
C ABOVE IS JUST A TRICK AS OTHERWISE GET 'UNBALANCED PARENTHESES' ERROR MESSAGE...
C
         VIXEN6=3.0D0*R(A,B)*R(B,C)*(R(A,B)**2*(Z(B)-Z(C)) + R(B,C)**2*(Z(B)-Z(A)))
         VIXEN7=U*(R(A,B)**2*(X(B)-X(C)) + R(B,C)**2*(X(B)-X(A)))
         VIXEN8= (VIXEN5*R(A,B)**3*R(B,C)**3 - VIXEN6*VIXEN7)/(R(A,B)**6*R(B,C)**6)
         VIXEN9=(2.0D0*X(B)-X(C)-X(A))*((Z(B)-Z(C))/R(B,C)**2 + (Z(B)-Z(A))/R(A,B)**2)/(R(A,B)*R(B,C))
         VIXEN4= VIXEN8 + VIXEN9

         ANGLEHELL(3*B-2,3*B)= ANGLEHELL(3*B-2,3*B) + VIXEN1*VIXEN4+VIXEN2*VIXEN3
         ANGLEHELL(3*B,3*B-2)= ANGLEHELL(3*B,3*B-2) + VIXEN1*VIXEN4+VIXEN2*VIXEN3
C
C DYBDZB,DZBDYB
C
         VIXEN2=U*(R(A,B)**2*(Y(B)-Y(C)) + R(B,C)**2*(Y(B)-Y(A)))/(R(A,B)**3*R(B,C)**3) - (2.0D0*Y(B)-Y(C)-Y(A))/(R(A,B)*R(B,C))
         VIXEN3=2.0D0*KTHETA*DTHETAZ2*(ST-(THETA(I)-THETA0)*CT)/ST**2

         VIXEN5=U*(2.0D0*(Z(B)-Z(A))*(Y(B)-Y(C))+2*(Z(B)-Z(C))*(Y(B)-Y(A)))
         VIXEN5= VIXEN5 + (2.0D0*Z(B)-Z(C)-Z(A))*(R(A,B)**2*(Y(B)-Y(C))+R(B,C)**2*(Y(B)-Y(A)))
C
C ABOVE IS JUST A TRICK AS OTHERWISE GET 'UNBALANCED PARENTHESES' ERROR MESSAGE...
C
         VIXEN6=3.0D0*R(A,B)*R(B,C)*(R(A,B)**2*(Z(B)-Z(C)) + R(B,C)**2*(Z(B)-Z(A)))
         VIXEN7=U*(R(A,B)**2*(Y(B)-Y(C)) + R(B,C)**2*(Y(B)-Y(A)))
         VIXEN8= (VIXEN5*R(A,B)**3*R(B,C)**3 - VIXEN6*VIXEN7)/(R(A,B)**6*R(B,C)**6)
         VIXEN9=(2.0D0*Y(B)-Y(C)-Y(A))*((Z(B)-Z(C))/R(B,C)**2 + (Z(B)-Z(A))/R(A,B)**2)/(R(A,B)*R(B,C))
         VIXEN4= VIXEN8 + VIXEN9

         ANGLEHELL(3*B-1,3*B)= ANGLEHELL(3*B-1,3*B) + VIXEN1*VIXEN4+VIXEN2*VIXEN3
         ANGLEHELL(3*B,3*B-1)= ANGLEHELL(3*B,3*B-1) + VIXEN1*VIXEN4+VIXEN2*VIXEN3

      END DO
C
C***************************
C TORSION ANGLE DERIVATIVES*
C***************************
C
      DO I=1,T
         A=DA1(I)
         B=DA2(I)
         C=DA3(I)
         D=DA4(I)
C         PRINT *,A,B,C,D,' : ',TYPECH(A),'-',TYPECH(B),'-',TYPECH(C),'-',TYPECH(D)

         RBC=R(B,C)
         RBC2=RBC**2
         RBC4=RBC**4
         RBC6=RBC**6
         XE=(X(C)-X(B))
         YE=(Y(C)-Y(B))
         ZE=(Z(C)-Z(B))
         LAMBDA=((X(B)-X(A))*XE+(Y(B)-Y(A))*YE+(Z(B)-Z(A))*ZE)/RBC2
         MU=    ((X(B)-X(D))*XE+(Y(B)-Y(D))*YE+(Z(B)-Z(D))*ZE)/RBC2
         XP= X(A)+LAMBDA*X(C)-(1+LAMBDA)*X(B)
         YP= Y(A)+LAMBDA*Y(C)-(1+LAMBDA)*Y(B)
         ZP= Z(A)+LAMBDA*Z(C)-(1+LAMBDA)*Z(B)
         XPP= X(D)+MU*X(C)-(1+MU)*X(B)
         YPP= Y(D)+MU*Y(C)-(1+MU)*Y(B)
         ZPP= Z(D)+MU*Z(C)-(1+MU)*Z(B)
         U= XP*XPP + YP*YPP + ZP*ZPP
         V= SQRT(XP**2 + YP**2 + ZP**2)*SQRT(XPP**2 + YPP**2 + ZPP**2)

         LAMBDACRAP= (X(B)-X(A))*XE + (Y(B)-Y(A))*YE + (Z(B)-Z(A))*ZE
         MUCRAP= (X(B)-X(D))*XE + (Y(B)-Y(D))*YE + (Z(B)-Z(D))*ZE

         DLAMBDAX1= (X(B)-X(C))/RBC2
         DLAMBDAY1= (Y(B)-Y(C))/RBC2
         DLAMBDAZ1= (Z(B)-Z(C))/RBC2
         DLAMBDAX2= (X(C)+X(A)-2.0D0*X(B))/RBC2 - 2.0D0*(X(B)-X(C))*LAMBDACRAP/RBC4 
         DLAMBDAY2= (Y(C)+Y(A)-2.0D0*Y(B))/RBC2 - 2.0D0*(Y(B)-Y(C))*LAMBDACRAP/RBC4
         DLAMBDAZ2= (Z(C)+Z(A)-2.0D0*Z(B))/RBC2 - 2.0D0*(Z(B)-Z(C))*LAMBDACRAP/RBC4
         DLAMBDAX3= (X(B)-X(A))/RBC2 - 2.0D0*XE*LAMBDACRAP/RBC4
         DLAMBDAY3= (Y(B)-Y(A))/RBC2 - 2.0D0*YE*LAMBDACRAP/RBC4
         DLAMBDAZ3= (Z(B)-Z(A))/RBC2 - 2.0D0*ZE*LAMBDACRAP/RBC4
         DMUX2= (X(C)+X(D)-2.0D0*X(B))/RBC2 - 2.0D0*(X(B)-X(C))*MUCRAP/RBC4
         DMUY2= (Y(C)+Y(D)-2.0D0*Y(B))/RBC2 - 2.0D0*(Y(B)-Y(C))*MUCRAP/RBC4
         DMUZ2= (Z(C)+Z(D)-2.0D0*Z(B))/RBC2 - 2.0D0*(Z(B)-Z(C))*MUCRAP/RBC4
         DMUX3= (X(B)-X(D))/RBC2 - 2.0D0*XE*MUCRAP/RBC4
         DMUY3= (Y(B)-Y(D))/RBC2 - 2.0D0*YE*MUCRAP/RBC4
         DMUZ3= (Z(B)-Z(D))/RBC2 - 2.0D0*ZE*MUCRAP/RBC4
         DMUX4= (X(B)-X(C))/RBC2
         DMUY4= (Y(B)-Y(C))/RBC2
         DMUZ4= (Z(B)-Z(C))/RBC2
C
C TESTING DERIVATIVE MATHS! FIRST THE DX1 TERMS
C
         VIXEN1=SQRT(XP**2 + YP**2 + ZP**2)
         VIXEN2=SQRT(XPP**2 + YPP**2 + ZPP**2)
         DPX1=(1.0D0/VIXEN1)*(XP*(1.0D0+XE*DLAMBDAX1) + YP*YE*DLAMBDAX1 + ZP*ZE*DLAMBDAX1)
         DVX1=VIXEN2*DPX1
         DUX1=XPP*(1.0D0+XE*DLAMBDAX1) + YPP*YE*DLAMBDAX1 + ZPP*ZE*DLAMBDAX1

         VIXEN5=(DVN(I)/DID(I))*4.0D0*COS(DPHI(I))*COS(DDELTA(I))
         TORSCRAP(3*A-2)= TORSCRAP(3*A-2)+ VIXEN5*(V*DUX1 - U*DVX1)/V**2

         DPY1=(1.0D0/VIXEN1)*(YP*(1.0D0+YE*DLAMBDAY1) + XP*XE*DLAMBDAY1 + ZP*ZE*DLAMBDAY1)
         DVY1=VIXEN2*DPY1
         DUY1=YPP*(1.0D0+YE*DLAMBDAY1) + XPP*XE*DLAMBDAY1 + ZPP*ZE*DLAMBDAY1

         TORSCRAP(3*A-1)= TORSCRAP(3*A-1)+ VIXEN5*(V*DUY1 - U*DVY1)/V**2

         DPZ1=(1.0D0/VIXEN1)*(ZP*(1.0D0+ZE*DLAMBDAZ1) + XP*XE*DLAMBDAZ1 + YP*YE*DLAMBDAZ1)
         DVZ1=VIXEN2*DPZ1
         DUZ1=ZPP*(1.0D0+ZE*DLAMBDAZ1) + XPP*XE*DLAMBDAZ1 + YPP*YE*DLAMBDAZ1

         TORSCRAP(3*A)= TORSCRAP(3*A)+ VIXEN5*(V*DUZ1 - U*DVZ1)/V**2
C
C NOW DX2 TERMS
C
         DUX2=XP*(XE*DMUX2-1.0D0-MU) + XPP*(XE*DLAMBDAX2-1.0D0-LAMBDA) + YE*(YP*DMUX2 + YPP*DLAMBDAX2) +
     1          ZE*(ZP*DMUX2 + ZPP*DLAMBDAX2)
         DPX2=(1.0D0/VIXEN1)*(XP*(XE*DLAMBDAX2-1.0D0-LAMBDA) + YP*YE*DLAMBDAX2 + ZP*ZE*DLAMBDAX2)
         VIXEN6=VIXEN2*DPX2
         DQX2=(1.0D0/VIXEN2)*(XPP*(XE*DMUX2-1.0D0-MU) + YPP*YE*DMUX2 + ZPP*ZE*DMUX2)
         VIXEN7=VIXEN1*DQX2
         DVX2=VIXEN6+VIXEN7

         TORSCRAP(3*B-2)= TORSCRAP(3*B-2)+ VIXEN5*(V*DUX2 - U*DVX2)/V**2

         DUY2=YP*(YE*DMUY2-1.0D0-MU) + YPP*(YE*DLAMBDAY2-1.0D0-LAMBDA) + XE*(XP*DMUY2 + XPP*DLAMBDAY2) +
     1          ZE*(ZP*DMUY2 + ZPP*DLAMBDAY2)
         DPY2=(1.0D0/VIXEN1)*(YP*(YE*DLAMBDAY2-1.0D0-LAMBDA) + XP*XE*DLAMBDAY2 + ZP*ZE*DLAMBDAY2)
         VIXEN6=VIXEN2*DPY2
         DQY2=(1.0D0/VIXEN2)*(YPP*(YE*DMUY2-1.0D0-MU) + XPP*XE*DMUY2 + ZPP*ZE*DMUY2)
         VIXEN7=VIXEN1*DQY2
         DVY2=VIXEN6+VIXEN7

         TORSCRAP(3*B-1)= TORSCRAP(3*B-1)+ VIXEN5*(V*DUY2 - U*DVY2)/V**2

         DUZ2=ZP*(ZE*DMUZ2-1.0D0-MU) + ZPP*(ZE*DLAMBDAZ2-1.0D0-LAMBDA) + XE*(XP*DMUZ2 + XPP*DLAMBDAZ2) +
     1          YE*(YP*DMUZ2 + YPP*DLAMBDAZ2)
         DPZ2=(1.0D0/VIXEN1)*(ZP*(ZE*DLAMBDAZ2-1.0D0-LAMBDA) + XP*XE*DLAMBDAZ2 + YP*YE*DLAMBDAZ2)
         VIXEN6=VIXEN2*DPZ2
         DQZ2=(1.0D0/VIXEN2)*(ZPP*(ZE*DMUZ2-1.0D0-MU) + XPP*XE*DMUZ2 + YPP*YE*DMUZ2)
         VIXEN7=VIXEN1*DQZ2
         DVZ2=VIXEN6+VIXEN7

         TORSCRAP(3*B)= TORSCRAP(3*B)+ VIXEN5*(V*DUZ2 - U*DVZ2)/V**2
C
C NOW DX3 TERMS
C
         DUX3=XP*(MU + XE*DMUX3) + XPP*(LAMBDA + XE*DLAMBDAX3) + YE*(YP*DMUX3 + YPP*DLAMBDAX3) +
     1          ZE*(ZP*DMUX3 + ZPP*DLAMBDAX3)
         DPX3=(1.0D0/VIXEN1)*(XP*(LAMBDA + XE*DLAMBDAX3) + YP*YE*DLAMBDAX3 + ZP*ZE*DLAMBDAX3)
         VIXEN6=VIXEN2*DPX3
         DQX3=(1.0D0/VIXEN2)*(XPP*(MU + XE*DMUX3) + YPP*YE*DMUX3 + ZPP*ZE*DMUX3)
         VIXEN7=VIXEN1*DQX3
         DVX3=VIXEN6+VIXEN7

         TORSCRAP(3*C-2)= TORSCRAP(3*C-2)+ VIXEN5*(V*DUX3 - U*DVX3)/V**2

         DUY3=YP*(MU + YE*DMUY3) + YPP*(LAMBDA + YE*DLAMBDAY3) + XE*(XP*DMUY3 + XPP*DLAMBDAY3) +
     1          ZE*(ZP*DMUY3 + ZPP*DLAMBDAY3)
         DPY3=(1.0D0/VIXEN1)*(YP*(LAMBDA + YE*DLAMBDAY3) + XP*XE*DLAMBDAY3 + ZP*ZE*DLAMBDAY3)
         VIXEN6=VIXEN2*DPY3
         DQY3=(1.0D0/VIXEN2)*(YPP*(MU + YE*DMUY3) + XPP*XE*DMUY3 + ZPP*ZE*DMUY3)
         VIXEN7=VIXEN1*DQY3
         DVY3=VIXEN6+VIXEN7

         TORSCRAP(3*C-1)= TORSCRAP(3*C-1)+ VIXEN5*(V*DUY3 - U*DVY3)/V**2

         DUZ3=ZP*(MU + ZE*DMUZ3) + ZPP*(LAMBDA + ZE*DLAMBDAZ3) + XE*(XP*DMUZ3 + XPP*DLAMBDAZ3) +
     1          YE*(YP*DMUZ3 + YPP*DLAMBDAZ3)
         DPZ3=(1.0D0/VIXEN1)*(ZP*(LAMBDA + ZE*DLAMBDAZ3) + XP*XE*DLAMBDAZ3 + YP*YE*DLAMBDAZ3)
         VIXEN6=VIXEN2*DPZ3
         DQZ3=(1.0D0/VIXEN2)*(ZPP*(MU + ZE*DMUZ3) + XPP*XE*DMUZ3 + YPP*YE*DMUZ3)
         VIXEN7=VIXEN1*DQZ3
         DVZ3=VIXEN6+VIXEN7

         TORSCRAP(3*C)= TORSCRAP(3*C)+ VIXEN5*(V*DUZ3 - U*DVZ3)/V**2
C
C NOW DX4 TERMS
C
         DUX4=XP*(1.0D0 + XE*DMUX4) + YP*YE*DMUX4 + ZP*ZE*DMUX4
         DQX4=(1.0D0/VIXEN2)*(XPP*(1.0D0 + XE*DMUX4) + YPP*YE*DMUX4 + ZPP*ZE*DMUX4)
         DVX4=VIXEN1*DQX4
  
         TORSCRAP(3*D-2)= TORSCRAP(3*D-2)+ VIXEN5*(V*DUX4 - U*DVX4)/V**2

         DUY4=YP*(1.0D0 + YE*DMUY4) + XP*XE*DMUY4 + ZP*ZE*DMUY4
         DQY4=(1.0D0/VIXEN2)*(YPP*(1.0D0 + YE*DMUY4) + XPP*XE*DMUY4 + ZPP*ZE*DMUY4)
         DVY4=VIXEN1*DQY4

         TORSCRAP(3*D-1)= TORSCRAP(3*D-1)+ VIXEN5*(V*DUY4 - U*DVY4)/V**2

         DUZ4=ZP*(1.0D0 + ZE*DMUZ4) + XP*XE*DMUZ4 + YP*YE*DMUZ4
         DQZ4=(1.0D0/VIXEN2)*(ZPP*(1.0D0 + ZE*DMUZ4) + XPP*XE*DMUZ4 + YPP*YE*DMUZ4)
         DVZ4=VIXEN1*DQZ4

         TORSCRAP(3*D)= TORSCRAP(3*D)+ VIXEN5*(V*DUZ4 - U*DVZ4)/V**2
C
C NOW SECOND DERIVATIVES
C
         CALL HAIRYHELL(U,V)
C
C DX1DX1
C
         N=(V*DUX1 - U*DVX1)/V**2
         VIXEN3= ((1.0D0+XE*DLAMBDAX1)**2 + YE**2*DLAMBDAX1**2 + ZE**2*DLAMBDAX1**2 - DPX1**2)/VIXEN1
         DNX1= -(U*VIXEN2/V**2)*VIXEN3 -2.0D0*DVX1*(V*DUX1 - U*DVX1)/V**3

         TORSHELL(3*A-2,3*A-2) =TORSHELL(3*A-2,3*A-2) + M*DNX1 + N*DMX1
C
C DY1DY1
C
         N=(V*DUY1 - U*DVY1)/V**2
         VIXEN3= ((1.0D0+YE*DLAMBDAY1)**2 + XE**2*DLAMBDAY1**2 + ZE**2*DLAMBDAY1**2 - DPY1**2)/VIXEN1
         DNY1= -(U*VIXEN2/V**2)*VIXEN3 -2.0D0*DVY1*(V*DUY1 - U*DVY1)/V**3

         TORSHELL(3*A-1,3*A-1) =TORSHELL(3*A-1,3*A-1) + M*DNY1 + N*DMY1
C
C DZ1DZ1
C
         N=(V*DUZ1 - U*DVZ1)/V**2
         VIXEN3= ((1.0D0+ZE*DLAMBDAZ1)**2 + XE**2*DLAMBDAZ1**2 + YE**2*DLAMBDAZ1**2 - DPZ1**2)/VIXEN1
         DNZ1= -(U*VIXEN2/V**2)*VIXEN3 -2.0D0*DVZ1*(V*DUZ1 - U*DVZ1)/V**3

         TORSHELL(3*A,3*A) =TORSHELL(3*A,3*A) + M*DNZ1 + N*DMZ1
C
C DX1DY1,DY1DX1
C
         N=(V*DUX1 - U*DVX1)/V**2
         VIXEN3= ((1.0D0+XE*DLAMBDAX1)*XE*DLAMBDAY1 + YE*DLAMBDAX1*(1.0D0+YE*DLAMBDAY1) 
     1          +ZE**2*DLAMBDAX1*DLAMBDAY1 - DPX1*DPY1)/VIXEN1
         DNY1= (DUX1*DVY1 - DUY1*DVX1 - U*VIXEN2*VIXEN3)/V**2 - 2.0D0*DVY1*(V*DUX1 - U*DVX1)/V**3
         
         TORSHELL(3*A-2,3*A-1)= TORSHELL(3*A-2,3*A-1) + M*DNY1 + N*DMY1
         TORSHELL(3*A-1,3*A-2)= TORSHELL(3*A-1,3*A-2) + M*DNY1 + N*DMY1
C
C DX1DZ1,DZ1DX1
C
         N=(V*DUX1 - U*DVX1)/V**2
         VIXEN3= ((1.0D0+XE*DLAMBDAX1)*XE*DLAMBDAZ1 + ZE*DLAMBDAX1*(1.0D0+ZE*DLAMBDAZ1)
     1          +YE**2*DLAMBDAX1*DLAMBDAZ1 - DPX1*DPZ1)/VIXEN1
         DNZ1= (DUX1*DVZ1 - DUZ1*DVX1 - U*VIXEN2*VIXEN3)/V**2 - 2.0D0*DVZ1*(V*DUX1 - U*DVX1)/V**3

         TORSHELL(3*A-2,3*A)= TORSHELL(3*A-2,3*A) + M*DNZ1 + N*DMZ1
         TORSHELL(3*A,3*A-2)= TORSHELL(3*A,3*A-2) + M*DNZ1 + N*DMZ1
C
C DY1DZ1,DZ1DY1
C
         N=(V*DUY1 - U*DVY1)/V**2
         VIXEN3= ((1.0D0+YE*DLAMBDAY1)*YE*DLAMBDAZ1 + ZE*DLAMBDAY1*(1.0D0+ZE*DLAMBDAZ1)
     1          +XE**2*DLAMBDAY1*DLAMBDAZ1 - DPY1*DPZ1)/VIXEN1
         DNZ1= (DUY1*DVZ1 - DUZ1*DVY1 - U*VIXEN2*VIXEN3)/V**2 - 2.0D0*DVZ1*(V*DUY1 - U*DVY1)/V**3

         TORSHELL(3*A-1,3*A)= TORSHELL(3*A-1,3*A) + M*DNZ1 + N*DMZ1
         TORSHELL(3*A,3*A-1)= TORSHELL(3*A,3*A-1) + M*DNZ1 + N*DMZ1
C
C DX1DX4,DX4DX1
C
         N=(V*DUX1 - U*DVX1)/V**2
         VIXEN3= (1.0D0+XE*DLAMBDAX1)*(1.0D0+XE*DMUX4) + (YE**2 + ZE**2)*DLAMBDAX1*DMUX4
         DNX4= (V*VIXEN3 + DVX4*DUX1 - U*DQX4*DPX1 - DVX1*DUX4)/V**2 - 2.0D0*DVX4*N/V

         TORSHELL(3*A-2,3*D-2)= TORSHELL(3*A-2,3*D-2) + M*DNX4 + N*DMX4
         TORSHELL(3*D-2,3*A-2)= TORSHELL(3*D-2,3*A-2) + M*DNX4 + N*DMX4
C
C DY1DY4,DY4DY1
C
         N=(V*DUY1 - U*DVY1)/V**2
         VIXEN3= (1.0D0+YE*DLAMBDAY1)*(1.0D0+YE*DMUY4) + (XE**2 + ZE**2)*DLAMBDAY1*DMUY4
         DNY4= (V*VIXEN3 + DVY4*DUY1 - U*DQY4*DPY1 - DVY1*DUY4)/V**2 - 2.0D0*DVY4*N/V

         TORSHELL(3*A-1,3*D-1)= TORSHELL(3*A-1,3*D-1) + M*DNY4 + N*DMY4
         TORSHELL(3*D-1,3*A-1)= TORSHELL(3*D-1,3*A-1) + M*DNY4 + N*DMY4
C
C DZ1DZ4,DZ4DZ1
C
         N=(V*DUZ1 - U*DVZ1)/V**2
         VIXEN3= (1.0D0+ZE*DLAMBDAZ1)*(1.0D0+ZE*DMUZ4) + (YE**2 + XE**2)*DLAMBDAZ1*DMUZ4
         DNZ4= (V*VIXEN3 + DVZ4*DUZ1 - U*DQZ4*DPZ1 - DVZ1*DUZ4)/V**2 - 2.0D0*DVZ4*N/V

         TORSHELL(3*A,3*D)= TORSHELL(3*A,3*D) + M*DNZ4 + N*DMZ4
         TORSHELL(3*D,3*A)= TORSHELL(3*D,3*A) + M*DNZ4 + N*DMZ4
C
C DX1DY4,DY4DX1
C
         N=(V*DUX1 - U*DVX1)/V**2
         VIXEN3= (1.0D0+XE*DLAMBDAX1)*XE*DMUY4 + YE*DLAMBDAX1*(1.0D0+YE*DMUY4) 
     1          + ZE**2*DLAMBDAX1*DMUY4
         DNY4= (V*VIXEN3 + DVY4*DUX1 - U*DQY4*DPX1 - DVX1*DUY4)/V**2 - 2.0D0*DVY4*N/V

         TORSHELL(3*A-2,3*D-1)= TORSHELL(3*A-2,3*D-1) + M*DNY4 + N*DMY4
         TORSHELL(3*D-1,3*A-2)= TORSHELL(3*D-1,3*A-2) + M*DNY4 + N*DMY4
C
C DX1DZ4,DZ4DX1
C
         N=(V*DUX1 - U*DVX1)/V**2
         VIXEN3= (1.0D0+XE*DLAMBDAX1)*XE*DMUZ4 + ZE*DLAMBDAX1*(1.0D0+ZE*DMUZ4)
     1          + YE**2*DLAMBDAX1*DMUZ4
         DNZ4= (V*VIXEN3 + DVZ4*DUX1 - U*DQZ4*DPX1 - DVX1*DUZ4)/V**2 - 2.0D0*DVZ4*N/V

         TORSHELL(3*A-2,3*D)= TORSHELL(3*A-2,3*D) + M*DNZ4 + N*DMZ4
         TORSHELL(3*D,3*A-2)= TORSHELL(3*D,3*A-2) + M*DNZ4 + N*DMZ4
C
C DY1DZ4,DZ4DY1
C
         N=(V*DUY1 - U*DVY1)/V**2
         VIXEN3= (1.0D0+YE*DLAMBDAY1)*YE*DMUZ4 + ZE*DLAMBDAY1*(1.0D0+ZE*DMUZ4)
     1          + XE**2*DLAMBDAY1*DMUZ4
         DNZ4= (V*VIXEN3 + DVZ4*DUY1 - U*DQZ4*DPY1 - DVY1*DUZ4)/V**2 - 2.0D0*DVZ4*N/V

         TORSHELL(3*A-1,3*D)= TORSHELL(3*A-1,3*D) + M*DNZ4 + N*DMZ4
         TORSHELL(3*D,3*A-1)= TORSHELL(3*D,3*A-1) + M*DNZ4 + N*DMZ4
C
C DY1DX4,DX4DY1
C
         N=(V*DUY1 - U*DVY1)/V**2
         VIXEN3= (1.0D0+YE*DLAMBDAY1)*YE*DMUX4 + XE*DLAMBDAY1*(1.0D0+XE*DMUX4)
     1          + ZE**2*DLAMBDAY1*DMUX4
         DNX4= (V*VIXEN3 + DVX4*DUY1 - U*DQX4*DPY1 - DVY1*DUX4)/V**2 - 2.0D0*DVX4*N/V

         TORSHELL(3*A-1,3*D-2)= TORSHELL(3*A-1,3*D-2) + M*DNX4 + N*DMX4
         TORSHELL(3*D-2,3*A-1)= TORSHELL(3*D-2,3*A-1) + M*DNX4 + N*DMX4
C
C DZ1DX4,DX4DZ1
C
         N=(V*DUZ1 - U*DVZ1)/V**2
         VIXEN3= (1.0D0+ZE*DLAMBDAZ1)*ZE*DMUX4 + XE*DLAMBDAZ1*(1.0D0+XE*DMUX4)
     1          + YE**2*DLAMBDAZ1*DMUX4
         DNX4= (V*VIXEN3 + DVX4*DUZ1 - U*DQX4*DPZ1 - DVZ1*DUX4)/V**2 - 2.0D0*DVX4*N/V

         TORSHELL(3*A,3*D-2)= TORSHELL(3*A,3*D-2) + M*DNX4 + N*DMX4
         TORSHELL(3*D-2,3*A)= TORSHELL(3*D-2,3*A) + M*DNX4 + N*DMX4
C
C DZ1DY4,DY4DZ1
C
         N=(V*DUZ1 - U*DVZ1)/V**2
         VIXEN3= (1.0D0+ZE*DLAMBDAZ1)*ZE*DMUY4 + YE*DLAMBDAZ1*(1.0D0+YE*DMUY4)
     1          + XE**2*DLAMBDAZ1*DMUY4
         DNY4= (V*VIXEN3 + DVY4*DUZ1 - U*DQY4*DPZ1 - DVZ1*DUY4)/V**2 - 2.0D0*DVY4*N/V

         TORSHELL(3*A,3*D-1)= TORSHELL(3*A,3*D-1) + M*DNY4 + N*DMY4
         TORSHELL(3*D-1,3*A)= TORSHELL(3*D-1,3*A) + M*DNY4 + N*DMY4
C
C DX4DX4
C
         N=(V*DUX4 - U*DVX4)/V**2
         VIXEN3= (VIXEN1/VIXEN2)*((1.0D0+XE*DMUX4)**2 + (YE**2 + ZE**2)*DMUX4**2 - DQX4**2)
         DNX4= -U*VIXEN3/V**2 - 2.0D0*DVX4*N/V

         TORSHELL(3*D-2,3*D-2)= TORSHELL(3*D-2,3*D-2) + M*DNX4 + N*DMX4
C
C DY4DY4
C
         N=(V*DUY4 - U*DVY4)/V**2
         VIXEN3= (VIXEN1/VIXEN2)*((1.0D0+YE*DMUY4)**2 + (XE**2 + ZE**2)*DMUY4**2 - DQY4**2)
         DNY4= -U*VIXEN3/V**2 - 2.0D0*DVY4*N/V

         TORSHELL(3*D-1,3*D-1)= TORSHELL(3*D-1,3*D-1) + M*DNY4 + N*DMY4
C
C DZ4DZ4
C
         N=(V*DUZ4 - U*DVZ4)/V**2
         VIXEN3= (VIXEN1/VIXEN2)*((1.0D0+ZE*DMUZ4)**2 + (YE**2 + XE**2)*DMUZ4**2 - DQZ4**2)
         DNZ4= -U*VIXEN3/V**2 - 2.0D0*DVZ4*N/V

         TORSHELL(3*D,3*D)= TORSHELL(3*D,3*D) + M*DNZ4 + N*DMZ4
C
C DX4DY4,DY4DX4
C
         N=(V*DUX4 - U*DVX4)/V**2
         VIXEN3= (VIXEN1/VIXEN2)*((1.0D0+XE*DMUX4)*XE*DMUY4 + YE*DMUX4*(1.0D0+YE*DMUY4)
     1          + ZE**2*DMUX4*DMUY4 - DQX4*DQY4)
         DNY4= (DVY4*DUX4 - U*VIXEN3 - DUY4*DVX4)/V**2 - 2.0D0*DVY4*N/V

         TORSHELL(3*D-2,3*D-1)= TORSHELL(3*D-2,3*D-1) + M*DNY4 + N*DMY4
         TORSHELL(3*D-1,3*D-2)= TORSHELL(3*D-1,3*D-2) + M*DNY4 + N*DMY4
C
C DX4DZ4,DZ4DX4
C
         N=(V*DUX4 - U*DVX4)/V**2
         VIXEN3= (VIXEN1/VIXEN2)*((1.0D0+XE*DMUX4)*XE*DMUZ4 + ZE*DMUX4*(1.0D0+ZE*DMUZ4)
     1          + YE**2*DMUX4*DMUZ4 - DQX4*DQZ4)
         DNZ4= (DVZ4*DUX4 - U*VIXEN3 - DUZ4*DVX4)/V**2 - 2.0D0*DVZ4*N/V

         TORSHELL(3*D-2,3*D)= TORSHELL(3*D-2,3*D) + M*DNZ4 + N*DMZ4
         TORSHELL(3*D,3*D-2)= TORSHELL(3*D,3*D-2) + M*DNZ4 + N*DMZ4
C
C DY4DZ4,DZ4DY4
C
         N=(V*DUY4 - U*DVY4)/V**2
         VIXEN3= (VIXEN1/VIXEN2)*((1.0D0+YE*DMUY4)*YE*DMUZ4 + ZE*DMUY4*(1.0D0+ZE*DMUZ4)
     1          + XE**2*DMUY4*DMUZ4 - DQY4*DQZ4)
         DNZ4= (DVZ4*DUY4 - U*VIXEN3 - DUZ4*DVY4)/V**2 - 2.0D0*DVZ4*N/V

         TORSHELL(3*D-1,3*D)= TORSHELL(3*D-1,3*D) + M*DNZ4 + N*DMZ4
         TORSHELL(3*D,3*D-1)= TORSHELL(3*D,3*D-1) + M*DNZ4 + N*DMZ4
C
C DX1DX2,DX2DX1
C
         N=(V*DUX1 - U*DVX1)/V**2
         VIXEN4= 1.0D0/RBC2 - 2.0D0*(X(B)-X(C))**2/RBC4
         VIXEN5= XP*(VIXEN4*XE-DLAMBDAX1) + (1.0D0+XE*DLAMBDAX1)*(XE*DLAMBDAX2-1.0D0-LAMBDA)
     1           + YP*YE*VIXEN4 + YE**2*DLAMBDAX1*DLAMBDAX2 + ZP*ZE*VIXEN4 
     2           + ZE**2*DLAMBDAX1*DLAMBDAX2 - DPX1*DPX2
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN5 + DPX1*DQX2 
         VIXEN6= XPP*(VIXEN4*XE-DLAMBDAX1) + (1.0D0+XE*DLAMBDAX1)*(XE*DMUX2-1.0D0-MU) +YPP*YE*VIXEN4
     1           + YE**2*DLAMBDAX1*DMUX2 + ZPP*ZE*VIXEN4 + ZE**2*DLAMBDAX1*DMUX2
         DNX2= (V*VIXEN6 + DVX2*DUX1 - U*VIXEN3 - DUX2*DVX1)/V**2 - 2.0D0*DVX2*N/V

         TORSHELL(3*A-2,3*B-2)= TORSHELL(3*A-2,3*B-2)+ M*DNX2 + N*DMX2
         TORSHELL(3*B-2,3*A-2)= TORSHELL(3*B-2,3*A-2)+ M*DNX2 + N*DMX2
C
C DY1DY2,DY2DY1
C
         N=(V*DUY1 - U*DVY1)/V**2
         VIXEN4= 1.0D0/RBC2 - 2.0D0*(Y(B)-Y(C))**2/RBC4
         VIXEN5= YP*(VIXEN4*YE-DLAMBDAY1) + (1.0D0+YE*DLAMBDAY1)*(YE*DLAMBDAY2-1.0D0-LAMBDA)
     1           + XP*XE*VIXEN4 + XE**2*DLAMBDAY1*DLAMBDAY2 + ZP*ZE*VIXEN4
     2           + ZE**2*DLAMBDAY1*DLAMBDAY2 - DPY1*DPY2
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN5 + DPY1*DQY2
         VIXEN6= YPP*(VIXEN4*YE-DLAMBDAY1) + (1.0D0+YE*DLAMBDAY1)*(YE*DMUY2-1.0D0-MU) +XPP*XE*VIXEN4
     1           + XE**2*DLAMBDAY1*DMUY2 + ZPP*ZE*VIXEN4 + ZE**2*DLAMBDAY1*DMUY2
         DNY2= (V*VIXEN6 + DVY2*DUY1 - U*VIXEN3 - DUY2*DVY1)/V**2 - 2.0D0*DVY2*N/V

         TORSHELL(3*A-1,3*B-1)= TORSHELL(3*A-1,3*B-1)+ M*DNY2 + N*DMY2
         TORSHELL(3*B-1,3*A-1)= TORSHELL(3*B-1,3*A-1)+ M*DNY2 + N*DMY2
C
C DZ1DZ2,DZ2DZ1
C
         N=(V*DUZ1 - U*DVZ1)/V**2
         VIXEN4= 1.0D0/RBC2 - 2.0D0*(Z(B)-Z(C))**2/RBC4
         VIXEN5= ZP*(VIXEN4*ZE-DLAMBDAZ1) + (1.0D0+ZE*DLAMBDAZ1)*(ZE*DLAMBDAZ2-1.0D0-LAMBDA)
     1           + XP*XE*VIXEN4 + XE**2*DLAMBDAZ1*DLAMBDAZ2 + YP*YE*VIXEN4
     2           + YE**2*DLAMBDAZ1*DLAMBDAZ2 - DPZ1*DPZ2
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN5 + DPZ1*DQZ2
         VIXEN6= ZPP*(VIXEN4*ZE-DLAMBDAZ1) + (1.0D0+ZE*DLAMBDAZ1)*(ZE*DMUZ2-1.0D0-MU) +XPP*XE*VIXEN4
     1           + XE**2*DLAMBDAZ1*DMUZ2 + YPP*YE*VIXEN4 + YE**2*DLAMBDAZ1*DMUZ2
         DNZ2= (V*VIXEN6 + DVZ2*DUZ1 - U*VIXEN3 - DUZ2*DVZ1)/V**2 - 2.0D0*DVZ2*N/V

         TORSHELL(3*A,3*B)= TORSHELL(3*A,3*B)+ M*DNZ2 + N*DMZ2
         TORSHELL(3*B,3*A)= TORSHELL(3*B,3*A)+ M*DNZ2 + N*DMZ2
C
C DX1DY2,DY2DX1
C
         N=(V*DUX1 - U*DVX1)/V**2
         VIXEN4= -2.0D0*(X(B)-X(C))*(Y(B)-Y(C))/RBC4
         VIXEN5= XP*XE*VIXEN4 + (1.0D0+XE*DLAMBDAX1)*XE*DLAMBDAY2 + YP*(YE*VIXEN4-DLAMBDAX1)
     1          + YE*DLAMBDAX1*(YE*DLAMBDAY2-1.0D0-LAMBDA) + ZP*ZE*VIXEN4 
     2          + ZE**2*DLAMBDAX1*DLAMBDAY2 - DPX1*DPY2
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN5 + DPX1*DQY2
         VIXEN6= XPP*XE*VIXEN4 + (1.0D0+XE*DLAMBDAX1)*XE*DMUY2 + YPP*(YE*VIXEN4-DLAMBDAX1)
     1          + YE*DLAMBDAX1*(YE*DMUY2-1.0D0-MU) + ZPP*ZE*VIXEN4 + ZE**2*DLAMBDAX1*DMUY2
         DNY2= (V*VIXEN6 + DVY2*DUX1 - U*VIXEN3 - DUY2*DVX1)/V**2 - 2.0D0*DVY2*N/V

         TORSHELL(3*A-2,3*B-1)= TORSHELL(3*A-2,3*B-1)+ M*DNY2 + N*DMY2
         TORSHELL(3*B-1,3*A-2)= TORSHELL(3*B-1,3*A-2)+ M*DNY2 + N*DMY2
C
C DX1DZ2,DZ2DX1
C
         N=(V*DUX1 - U*DVX1)/V**2
         VIXEN4= -2.0D0*(X(B)-X(C))*(Z(B)-Z(C))/RBC4
         VIXEN5= XP*XE*VIXEN4 + (1.0D0+XE*DLAMBDAX1)*XE*DLAMBDAZ2 + ZP*(ZE*VIXEN4-DLAMBDAX1)
     1          + ZE*DLAMBDAX1*(ZE*DLAMBDAZ2-1.0D0-LAMBDA) + YP*YE*VIXEN4
     2          + YE**2*DLAMBDAX1*DLAMBDAZ2 - DPX1*DPZ2
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN5 + DPX1*DQZ2
         VIXEN6= XPP*XE*VIXEN4 + (1.0D0+XE*DLAMBDAX1)*XE*DMUZ2 + ZPP*(ZE*VIXEN4-DLAMBDAX1)
     1          + ZE*DLAMBDAX1*(ZE*DMUZ2-1.0D0-MU) + YPP*YE*VIXEN4 + YE**2*DLAMBDAX1*DMUZ2
         DNZ2= (V*VIXEN6 + DVZ2*DUX1 - U*VIXEN3 - DUZ2*DVX1)/V**2 - 2.0D0*DVZ2*N/V

         TORSHELL(3*A-2,3*B)= TORSHELL(3*A-2,3*B)+ M*DNZ2 + N*DMZ2
         TORSHELL(3*B,3*A-2)= TORSHELL(3*B,3*A-2)+ M*DNZ2 + N*DMZ2
C
C DY1DZ2,DZ2DY1
C
         N=(V*DUY1 - U*DVY1)/V**2
         VIXEN4= -2.0D0*(Y(B)-Y(C))*(Z(B)-Z(C))/RBC4
         VIXEN5= YP*YE*VIXEN4 + (1.0D0+YE*DLAMBDAY1)*YE*DLAMBDAZ2 + ZP*(ZE*VIXEN4-DLAMBDAY1)
     1          + ZE*DLAMBDAY1*(ZE*DLAMBDAZ2-1.0D0-LAMBDA) + XP*XE*VIXEN4
     2          + XE**2*DLAMBDAY1*DLAMBDAZ2 - DPY1*DPZ2
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN5 + DPY1*DQZ2
         VIXEN6= YPP*YE*VIXEN4 + (1.0D0+YE*DLAMBDAY1)*YE*DMUZ2 + ZPP*(ZE*VIXEN4-DLAMBDAY1)
     1          + ZE*DLAMBDAY1*(ZE*DMUZ2-1.0D0-MU) + XPP*XE*VIXEN4 + XE**2*DLAMBDAY1*DMUZ2
         DNZ2= (V*VIXEN6 + DVZ2*DUY1 - U*VIXEN3 - DUZ2*DVY1)/V**2 - 2.0D0*DVZ2*N/V

         TORSHELL(3*A-1,3*B)= TORSHELL(3*A-1,3*B)+ M*DNZ2 + N*DMZ2
         TORSHELL(3*B,3*A-1)= TORSHELL(3*B,3*A-1)+ M*DNZ2 + N*DMZ2
C
C DY1DX2,DX2DY1
C
         N=(V*DUY1 - U*DVY1)/V**2
         VIXEN4= -2.0D0*(Y(B)-Y(C))*(X(B)-X(C))/RBC4
         VIXEN5= YP*YE*VIXEN4 + (1.0D0+YE*DLAMBDAY1)*YE*DLAMBDAX2 + XP*(XE*VIXEN4-DLAMBDAY1)
     1          + XE*DLAMBDAY1*(XE*DLAMBDAX2-1.0D0-LAMBDA) + ZP*ZE*VIXEN4
     2          + ZE**2*DLAMBDAY1*DLAMBDAX2 - DPY1*DPX2
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN5 + DPY1*DQX2
         VIXEN6= YPP*YE*VIXEN4 + (1.0D0+YE*DLAMBDAY1)*YE*DMUX2 + XPP*(XE*VIXEN4-DLAMBDAY1)
     1          + XE*DLAMBDAY1*(XE*DMUX2-1.0D0-MU) + ZPP*ZE*VIXEN4 + ZE**2*DLAMBDAY1*DMUX2
         DNX2= (V*VIXEN6 + DVX2*DUY1 - U*VIXEN3 - DUX2*DVY1)/V**2 - 2.0D0*DVX2*N/V

         TORSHELL(3*A-1,3*B-2)= TORSHELL(3*A-1,3*B-2)+ M*DNX2 + N*DMX2
         TORSHELL(3*B-2,3*A-1)= TORSHELL(3*B-2,3*A-1)+ M*DNX2 + N*DMX2
C
C DZ1DX2,DX2DZ1
C
         N=(V*DUZ1 - U*DVZ1)/V**2
         VIXEN4= -2.0D0*(Z(B)-Z(C))*(X(B)-X(C))/RBC4
         VIXEN5= ZP*ZE*VIXEN4 + (1.0D0+ZE*DLAMBDAZ1)*ZE*DLAMBDAX2 + XP*(XE*VIXEN4-DLAMBDAZ1)
     1          + XE*DLAMBDAZ1*(XE*DLAMBDAX2-1.0D0-LAMBDA) + YP*YE*VIXEN4
     2          + YE**2*DLAMBDAZ1*DLAMBDAX2 - DPZ1*DPX2
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN5 + DPZ1*DQX2
         VIXEN6= ZPP*ZE*VIXEN4 + (1.0D0+ZE*DLAMBDAZ1)*ZE*DMUX2 + XPP*(XE*VIXEN4-DLAMBDAZ1)
     1          + XE*DLAMBDAZ1*(XE*DMUX2-1-MU) + YPP*YE*VIXEN4 + YE**2*DLAMBDAZ1*DMUX2
         DNX2= (V*VIXEN6 + DVX2*DUZ1 - U*VIXEN3 - DUX2*DVZ1)/V**2 - 2.0D0*DVX2*N/V

         TORSHELL(3*A,3*B-2)= TORSHELL(3*A,3*B-2)+ M*DNX2 + N*DMX2
         TORSHELL(3*B-2,3*A)= TORSHELL(3*B-2,3*A)+ M*DNX2 + N*DMX2
C
C DZ1DY2,DY2DZ1
C
         N=(V*DUZ1 - U*DVZ1)/V**2
         VIXEN4= -2.0D0*(Z(B)-Z(C))*(Y(B)-Y(C))/RBC4
         VIXEN5= ZP*ZE*VIXEN4 + (1.0D0+ZE*DLAMBDAZ1)*ZE*DLAMBDAY2 + YP*(YE*VIXEN4-DLAMBDAZ1)
     1          + YE*DLAMBDAZ1*(YE*DLAMBDAY2-1.0D0-LAMBDA) + XP*XE*VIXEN4
     2          + XE**2*DLAMBDAZ1*DLAMBDAY2 - DPZ1*DPY2
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN5 + DPZ1*DQY2
         VIXEN6= ZPP*ZE*VIXEN4 + (1.0D0+ZE*DLAMBDAZ1)*ZE*DMUY2 + YPP*(YE*VIXEN4-DLAMBDAZ1)
     1          + YE*DLAMBDAZ1*(YE*DMUY2-1-MU) + XPP*XE*VIXEN4 + XE**2*DLAMBDAZ1*DMUY2
         DNY2= (V*VIXEN6 + DVY2*DUZ1 - U*VIXEN3 - DUY2*DVZ1)/V**2 - 2.0D0*DVY2*N/V

         TORSHELL(3*A,3*B-1)= TORSHELL(3*A,3*B-1)+ M*DNY2 + N*DMY2
         TORSHELL(3*B-1,3*A)= TORSHELL(3*B-1,3*A)+ M*DNY2 + N*DMY2
C
C DX1DX3,DX3DX1
C
         N=(V*DUX1 - U*DVX1)/V**2
         VIXEN4= 2.0D0*(X(B)-X(C))**2/RBC4 - 1/RBC2
         VIXEN5= XP*(XE*VIXEN4+DLAMBDAX1) + (1.0D0+XE*DLAMBDAX1)*(XE*DLAMBDAX3+LAMBDA)
     1          + (YP*YE+ZP*ZE)*VIXEN4 + (YE**2+ZE**2)*DLAMBDAX1*DLAMBDAX3 - DPX1*DPX3
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN5 + DPX1*DQX3
         VIXEN6= XPP*(XE*VIXEN4+DLAMBDAX1) + (1.0D0+XE*DLAMBDAX1)*(XE*DMUX3+MU)
     1          + (YPP*YE+ZPP*ZE)*VIXEN4 + (YE**2+ZE**2)*DLAMBDAX1*DMUX3
         DNX3= (V*VIXEN6 + DVX3*DUX1 - U*VIXEN3 - DUX3*DVX1)/V**2 - 2.0D0*DVX3*N/V

         TORSHELL(3*A-2,3*C-2)= TORSHELL(3*A-2,3*C-2)+ M*DNX3 + N*DMX3
         TORSHELL(3*C-2,3*A-2)= TORSHELL(3*C-2,3*A-2)+ M*DNX3 + N*DMX3
C
C DY1DY3,DY3DY1
C
         N=(V*DUY1 - U*DVY1)/V**2
         VIXEN4= 2.0D0*(Y(B)-Y(C))**2/RBC4 - 1/RBC2
         VIXEN5= YP*(YE*VIXEN4+DLAMBDAY1) + (1.0D0+YE*DLAMBDAY1)*(YE*DLAMBDAY3+LAMBDA)
     1          + (XP*XE+ZP*ZE)*VIXEN4 + (XE**2+ZE**2)*DLAMBDAY1*DLAMBDAY3 - DPY1*DPY3
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN5 + DPY1*DQY3
         VIXEN6= YPP*(YE*VIXEN4+DLAMBDAY1) + (1.0D0+YE*DLAMBDAY1)*(YE*DMUY3+MU)
     1          + (XPP*XE+ZPP*ZE)*VIXEN4 + (XE**2+ZE**2)*DLAMBDAY1*DMUY3
         DNY3= (V*VIXEN6 + DVY3*DUY1 - U*VIXEN3 - DUY3*DVY1)/V**2 - 2.0D0*DVY3*N/V

         TORSHELL(3*A-1,3*C-1)= TORSHELL(3*A-1,3*C-1)+ M*DNY3 + N*DMY3
         TORSHELL(3*C-1,3*A-1)= TORSHELL(3*C-1,3*A-1)+ M*DNY3 + N*DMY3
C
C DZ1DZ3,DZ3DZ1
C
         N=(V*DUZ1 - U*DVZ1)/V**2
         VIXEN4= 2.0D0*(Z(B)-Z(C))**2/RBC4 - 1/RBC2
         VIXEN5= ZP*(ZE*VIXEN4+DLAMBDAZ1) + (1.0D0+ZE*DLAMBDAZ1)*(ZE*DLAMBDAZ3+LAMBDA)
     1          + (XP*XE+YP*YE)*VIXEN4 + (XE**2+YE**2)*DLAMBDAZ1*DLAMBDAZ3 - DPZ1*DPZ3
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN5 + DPZ1*DQZ3
         VIXEN6= ZPP*(ZE*VIXEN4+DLAMBDAZ1) + (1.0D0+ZE*DLAMBDAZ1)*(ZE*DMUZ3+MU)
     1          + (XPP*XE+YPP*YE)*VIXEN4 + (XE**2+YE**2)*DLAMBDAZ1*DMUZ3
         DNZ3= (V*VIXEN6 + DVZ3*DUZ1 - U*VIXEN3 - DUZ3*DVZ1)/V**2 - 2.0D0*DVZ3*N/V

         TORSHELL(3*A,3*C)= TORSHELL(3*A,3*C)+ M*DNZ3 + N*DMZ3
         TORSHELL(3*C,3*A)= TORSHELL(3*C,3*A)+ M*DNZ3 + N*DMZ3
C
C DX1DY3
C
         N=(V*DUX1 - U*DVX1)/V**2
         VIXEN4= 2.0D0*(X(B)-X(C))*(Y(B)-Y(C))/RBC4
         VIXEN5= XP*XE*VIXEN4 + (1.0D0+XE*DLAMBDAX1)*XE*DLAMBDAY3 
     1          + YP*(YE*VIXEN4+DLAMBDAX1) + YE*DLAMBDAX1*(YE*DLAMBDAY3+LAMBDA) 
     2          + ZP*ZE*VIXEN4 + ZE**2*DLAMBDAX1*DLAMBDAY3 - DPX1*DPY3
         VIXEN3=(VIXEN2/VIXEN1)*VIXEN5 + DPX1*DQY3
         VIXEN6= XPP*XE*VIXEN4 + (1.0D0+XE*DLAMBDAX1)*XE*DMUY3 + YPP*(YE*VIXEN4+DLAMBDAX1)
     1          + YE*DLAMBDAX1*(YE*DMUY3+MU) + ZPP*ZE*VIXEN4 + ZE**2*DLAMBDAX1*DMUY3
         DNY3= (V*VIXEN6 + DVY3*DUX1 - U*VIXEN3 - DUY3*DVX1)/V**2 - 2.0D0*DVY3*N/V

         TORSHELL(3*A-2,3*C-1)= TORSHELL(3*A-2,3*C-1)+ M*DNY3 + N*DMY3
         TORSHELL(3*C-1,3*A-2)= TORSHELL(3*C-1,3*A-2)+ M*DNY3 + N*DMY3
C
C DX1DZ3
C
         N=(V*DUX1 - U*DVX1)/V**2
         VIXEN4= 2.0D0*(X(B)-X(C))*(Z(B)-Z(C))/RBC4
         VIXEN5= XP*XE*VIXEN4 + (1.0D0+XE*DLAMBDAX1)*XE*DLAMBDAZ3
     1          + ZP*(ZE*VIXEN4+DLAMBDAX1) + ZE*DLAMBDAX1*(ZE*DLAMBDAZ3+LAMBDA)
     2          + YP*YE*VIXEN4 + YE**2*DLAMBDAX1*DLAMBDAZ3 - DPX1*DPZ3
         VIXEN3=(VIXEN2/VIXEN1)*VIXEN5 + DPX1*DQZ3
         VIXEN6= XPP*XE*VIXEN4 + (1.0D0+XE*DLAMBDAX1)*XE*DMUZ3 + ZPP*(ZE*VIXEN4+DLAMBDAX1)
     1          + ZE*DLAMBDAX1*(ZE*DMUZ3+MU) + YPP*YE*VIXEN4 + YE**2*DLAMBDAX1*DMUZ3
         DNZ3= (V*VIXEN6 + DVZ3*DUX1 - U*VIXEN3 - DUZ3*DVX1)/V**2 - 2.0D0*DVZ3*N/V

         TORSHELL(3*A-2,3*C)= TORSHELL(3*A-2,3*C)+ M*DNZ3 + N*DMZ3
         TORSHELL(3*C,3*A-2)= TORSHELL(3*C,3*A-2)+ M*DNZ3 + N*DMZ3
C
C DY1DZ3
C
         N=(V*DUY1 - U*DVY1)/V**2
         VIXEN4= 2.0D0*(Y(B)-Y(C))*(Z(B)-Z(C))/RBC4
         VIXEN5= YP*YE*VIXEN4 + (1.0D0+YE*DLAMBDAY1)*YE*DLAMBDAZ3
     1          + ZP*(ZE*VIXEN4+DLAMBDAY1) + ZE*DLAMBDAY1*(ZE*DLAMBDAZ3+LAMBDA)
     2          + XP*XE*VIXEN4 + XE**2*DLAMBDAY1*DLAMBDAZ3 - DPY1*DPZ3
         VIXEN3=(VIXEN2/VIXEN1)*VIXEN5 + DPY1*DQZ3
         VIXEN6= YPP*YE*VIXEN4 + (1.0D0+YE*DLAMBDAY1)*YE*DMUZ3 + ZPP*(ZE*VIXEN4+DLAMBDAY1)
     1          + ZE*DLAMBDAY1*(ZE*DMUZ3+MU) + XPP*XE*VIXEN4 + XE**2*DLAMBDAY1*DMUZ3
         DNZ3= (V*VIXEN6 + DVZ3*DUY1 - U*VIXEN3 - DUZ3*DVY1)/V**2 - 2.0D0*DVZ3*N/V

         TORSHELL(3*A-1,3*C)= TORSHELL(3*A-1,3*C)+ M*DNZ3 + N*DMZ3
         TORSHELL(3*C,3*A-1)= TORSHELL(3*C,3*A-1)+ M*DNZ3 + N*DMZ3
C
C DY1DX3
C
         N=(V*DUY1 - U*DVY1)/V**2
         VIXEN4= 2.0D0*(Y(B)-Y(C))*(X(B)-X(C))/RBC4
         VIXEN5= YP*YE*VIXEN4 + (1.0D0+YE*DLAMBDAY1)*YE*DLAMBDAX3
     1          + XP*(XE*VIXEN4+DLAMBDAY1) + XE*DLAMBDAY1*(XE*DLAMBDAX3+LAMBDA)
     2          + ZP*ZE*VIXEN4 + ZE**2*DLAMBDAY1*DLAMBDAX3 - DPY1*DPX3
         VIXEN3=(VIXEN2/VIXEN1)*VIXEN5 + DPY1*DQX3
         VIXEN6= YPP*YE*VIXEN4 + (1.0D0+YE*DLAMBDAY1)*YE*DMUX3 + XPP*(XE*VIXEN4+DLAMBDAY1)
     1          + XE*DLAMBDAY1*(XE*DMUX3+MU) + ZPP*ZE*VIXEN4 + ZE**2*DLAMBDAY1*DMUX3
         DNX3= (V*VIXEN6 + DVX3*DUY1 - U*VIXEN3 - DUX3*DVY1)/V**2 - 2.0D0*DVX3*N/V

         TORSHELL(3*A-1,3*C-2)= TORSHELL(3*A-1,3*C-2)+ M*DNX3 + N*DMX3
         TORSHELL(3*C-2,3*A-1)= TORSHELL(3*C-2,3*A-1)+ M*DNX3 + N*DMX3
C
C DZ1DX3
C
         N=(V*DUZ1 - U*DVZ1)/V**2
         VIXEN4= 2.0D0*(Z(B)-Z(C))*(X(B)-X(C))/RBC4
         VIXEN5= ZP*ZE*VIXEN4 + (1.0D0+ZE*DLAMBDAZ1)*ZE*DLAMBDAX3
     1          + XP*(XE*VIXEN4+DLAMBDAZ1) + XE*DLAMBDAZ1*(XE*DLAMBDAX3+LAMBDA)
     2          + YP*YE*VIXEN4 + YE**2*DLAMBDAZ1*DLAMBDAX3 - DPZ1*DPX3
         VIXEN3=(VIXEN2/VIXEN1)*VIXEN5 + DPZ1*DQX3
         VIXEN6= ZPP*ZE*VIXEN4 + (1.0D0+ZE*DLAMBDAZ1)*ZE*DMUX3 + XPP*(XE*VIXEN4+DLAMBDAZ1)
     1          + XE*DLAMBDAZ1*(XE*DMUX3+MU) + YPP*YE*VIXEN4 + YE**2*DLAMBDAZ1*DMUX3
         DNX3= (V*VIXEN6 + DVX3*DUZ1 - U*VIXEN3 - DUX3*DVZ1)/V**2 - 2.0D0*DVX3*N/V

         TORSHELL(3*A,3*C-2)= TORSHELL(3*A,3*C-2)+ M*DNX3 + N*DMX3
         TORSHELL(3*C-2,3*A)= TORSHELL(3*C-2,3*A)+ M*DNX3 + N*DMX3
C
C DZ1DY3
C
         N=(V*DUZ1 - U*DVZ1)/V**2
         VIXEN4= 2.0D0*(Z(B)-Z(C))*(Y(B)-Y(C))/RBC4
         VIXEN5= ZP*ZE*VIXEN4 + (1.0D0+ZE*DLAMBDAZ1)*ZE*DLAMBDAY3
     1          + YP*(YE*VIXEN4+DLAMBDAZ1) + YE*DLAMBDAZ1*(YE*DLAMBDAY3+LAMBDA)
     2          + XP*XE*VIXEN4 + XE**2*DLAMBDAZ1*DLAMBDAY3 - DPZ1*DPY3
         VIXEN3=(VIXEN2/VIXEN1)*VIXEN5 + DPZ1*DQY3
         VIXEN6= ZPP*ZE*VIXEN4 + (1.0D0+ZE*DLAMBDAZ1)*ZE*DMUY3 + YPP*(YE*VIXEN4+DLAMBDAZ1)
     1          + YE*DLAMBDAZ1*(YE*DMUY3+MU) + XPP*XE*VIXEN4 + XE**2*DLAMBDAZ1*DMUY3
         DNY3= (V*VIXEN6 + DVY3*DUZ1 - U*VIXEN3 - DUY3*DVZ1)/V**2 - 2.0D0*DVY3*N/V

         TORSHELL(3*A,3*C-1)= TORSHELL(3*A,3*C-1)+ M*DNY3 + N*DMY3
         TORSHELL(3*C-1,3*A)= TORSHELL(3*C-1,3*A)+ M*DNY3 + N*DMY3
C
C DX4DX2,DX2DX4
C
         N=(V*DUX4 - U*DVX4)/V**2
         VIXEN4= 1.0D0/RBC2 - 2.0D0*(X(B)-X(C))**2/RBC4
         VIXEN5= XPP*(XE*VIXEN4-DMUX4) + (1.0D0+XE*DMUX4)
     1          *(XE*DMUX2-1.0D0-MU) + YPP*YE*VIXEN4
     2          + YE**2*DMUX4*DMUX2 + ZPP*ZE*VIXEN4
     3          + ZE**2*DMUX4*DMUX2 - DQX4*DQX2
         VIXEN3=(VIXEN1/VIXEN2)*VIXEN5 + DPX2*DQX4
         VIXEN6= XP*(XE*VIXEN4-DMUX4) + (1.0D0+XE*DMUX4)
     1          *(XE*DLAMBDAX2-1.0D0-LAMBDA) + YP*YE*VIXEN4
     2          + YE**2*DMUX4*DLAMBDAX2 + ZP*ZE*VIXEN4
     3          + ZE**2*DMUX4*DLAMBDAX2
         DNX2= (V*VIXEN6 + DVX2*DUX4 - U*VIXEN3 - DUX2*DVX4)/V**2
     1        - 2.0D0*DVX2*(V*DUX4 - U*DVX4)/V**3

         TORSHELL(3*D-2,3*B-2)= TORSHELL(3*D-2,3*B-2)+ M*DNX2 + N*DMX2
         TORSHELL(3*B-2,3*D-2)= TORSHELL(3*B-2,3*D-2)+ M*DNX2 + N*DMX2
C
C DY4DY2,DY2DY4
C
         N=(V*DUY4 - U*DVY4)/V**2
         VIXEN4= 1.0D0/RBC2 - 2.0D0*(Y(B)-Y(C))**2/RBC4
         VIXEN5= YPP*(YE*VIXEN4-DMUY4) + (1.0D0+YE*DMUY4)
     1          *(YE*DMUY2-1.0D0-MU) + XPP*XE*VIXEN4
     2          + XE**2*DMUY4*DMUY2 + ZPP*ZE*VIXEN4
     3          + ZE**2*DMUY4*DMUY2 - DQY4*DQY2
         VIXEN3=(VIXEN1/VIXEN2)*VIXEN5 + DPY2*DQY4
         VIXEN6= YP*(YE*VIXEN4-DMUY4) + (1.0D0+YE*DMUY4)
     1          *(YE*DLAMBDAY2-1.0D0-LAMBDA) + XP*XE*VIXEN4
     2          + XE**2*DMUY4*DLAMBDAY2 + ZP*ZE*VIXEN4
     3          + ZE**2*DMUY4*DLAMBDAY2
         DNY2= (V*VIXEN6 + DVY2*DUY4 - U*VIXEN3 - DUY2*DVY4)/V**2
     1        - 2.0D0*DVY2*(V*DUY4 - U*DVY4)/V**3

         TORSHELL(3*D-1,3*B-1)= TORSHELL(3*D-1,3*B-1)+ M*DNY2 + N*DMY2
         TORSHELL(3*B-1,3*D-1)= TORSHELL(3*B-1,3*D-1)+ M*DNY2 + N*DMY2
C
C DZ4DZ2,DZ2DZ4
C
         N=(V*DUZ4 - U*DVZ4)/V**2
         VIXEN4= 1.0D0/RBC2 - 2.0D0*(Z(B)-Z(C))**2/RBC4
         VIXEN5= ZPP*(ZE*VIXEN4-DMUZ4) + (1.0D0+ZE*DMUZ4)
     1          *(ZE*DMUZ2-1.0D0-MU) + XPP*XE*VIXEN4
     2          + XE**2*DMUZ4*DMUZ2 + YPP*YE*VIXEN4
     3          + YE**2*DMUZ4*DMUZ2 - DQZ4*DQZ2
         VIXEN3=(VIXEN1/VIXEN2)*VIXEN5 + DPZ2*DQZ4
         VIXEN6= ZP*(ZE*VIXEN4-DMUZ4) + (1.0D0+ZE*DMUZ4)
     1          *(ZE*DLAMBDAZ2-1.0D0-LAMBDA) + XP*XE*VIXEN4
     2          + XE**2*DMUZ4*DLAMBDAZ2 + YP*YE*VIXEN4
     3          + YE**2*DMUZ4*DLAMBDAZ2
         DNZ2= (V*VIXEN6 + DVZ2*DUZ4 - U*VIXEN3 - DUZ2*DVZ4)/V**2
     1        - 2.0D0*DVZ2*(V*DUZ4 - U*DVZ4)/V**3

         TORSHELL(3*D,3*B)= TORSHELL(3*D,3*B)+ M*DNZ2 + N*DMZ2
         TORSHELL(3*B,3*D)= TORSHELL(3*B,3*D)+ M*DNZ2 + N*DMZ2
C
C DX4DY2,DY2DX4
C
         N=(V*DUX4 - U*DVX4)/V**2
         VIXEN4= -2.0D0*(X(B)-X(C))*(Y(B)-Y(C))/RBC4
         VIXEN5= XPP*XE*VIXEN4 + (1.0D0+XE*DMUX4)*XE*DMUY2 + YPP*(YE*VIXEN4-DMUX4)
     1          + YE*DMUX4*(YE*DMUY2-1.0D0-MU) + ZPP*ZE*VIXEN4 + ZE**2*DMUX4*DMUY2 - DQX4*DQY2
         VIXEN3=(VIXEN1/VIXEN2)*VIXEN5 + DPY2*DQX4
         VIXEN6= XP*XE*VIXEN4 + (1.0D0+XE*DMUX4)*XE*DLAMBDAY2 + YP*(YE*VIXEN4-DMUX4)
     1          + YE*DMUX4*(YE*DLAMBDAY2-1.0D0-LAMBDA) + ZP*ZE*VIXEN4 + ZE**2*DMUX4*DLAMBDAY2
         DNY2= (V*VIXEN6 + DVY2*DUX4 - U*VIXEN3 - DUY2*DVX4)/V**2 - 2.0D0*DVY2*N/V

         TORSHELL(3*D-2,3*B-1)= TORSHELL(3*D-2,3*B-1)+ M*DNY2 + N*DMY2
         TORSHELL(3*B-1,3*D-2)= TORSHELL(3*B-1,3*D-2)+ M*DNY2 + N*DMY2
C
C DX4DZ2,DZ2DX4
C
         N=(V*DUX4 - U*DVX4)/V**2
         VIXEN4= -2.0D0*(X(B)-X(C))*(Z(B)-Z(C))/RBC4
         VIXEN5= XPP*XE*VIXEN4 + (1.0D0+XE*DMUX4)*XE*DMUZ2 + ZPP*(ZE*VIXEN4-DMUX4)
     1          + ZE*DMUX4*(ZE*DMUZ2-1.0D0-MU) + YPP*YE*VIXEN4 + YE**2*DMUX4*DMUZ2 - DQX4*DQZ2
         VIXEN3=(VIXEN1/VIXEN2)*VIXEN5 + DPZ2*DQX4
         VIXEN6= XP*XE*VIXEN4 + (1.0D0+XE*DMUX4)*XE*DLAMBDAZ2 + ZP*(ZE*VIXEN4-DMUX4)
     1          + ZE*DMUX4*(ZE*DLAMBDAZ2-1.0D0-LAMBDA) + YP*YE*VIXEN4 + YE**2*DMUX4*DLAMBDAZ2
         DNZ2= (V*VIXEN6 + DVZ2*DUX4 - U*VIXEN3 - DUZ2*DVX4)/V**2 - 2.0D0*DVZ2*N/V

         TORSHELL(3*D-2,3*B)= TORSHELL(3*D-2,3*B)+ M*DNZ2 + N*DMZ2
         TORSHELL(3*B,3*D-2)= TORSHELL(3*B,3*D-2)+ M*DNZ2 + N*DMZ2
C
C DY4DZ2,DZ2DY4
C
         N=(V*DUY4 - U*DVY4)/V**2
         VIXEN4= -2.0D0*(Y(B)-Y(C))*(Z(B)-Z(C))/RBC4
         VIXEN5= YPP*YE*VIXEN4 + (1.0D0+YE*DMUY4)*YE*DMUZ2 + ZPP*(ZE*VIXEN4-DMUY4)
     1          + ZE*DMUY4*(ZE*DMUZ2-1.0D0-MU) + XPP*XE*VIXEN4 + XE**2*DMUY4*DMUZ2 - DQY4*DQZ2
         VIXEN3=(VIXEN1/VIXEN2)*VIXEN5 + DPZ2*DQY4
         VIXEN6= YP*YE*VIXEN4 + (1.0D0+YE*DMUY4)*YE*DLAMBDAZ2 + ZP*(ZE*VIXEN4-DMUY4)
     1          + ZE*DMUY4*(ZE*DLAMBDAZ2-1.0D0-LAMBDA) + XP*XE*VIXEN4 + XE**2*DMUY4*DLAMBDAZ2
         DNZ2= (V*VIXEN6 + DVZ2*DUY4 - U*VIXEN3 - DUZ2*DVY4)/V**2 - 2.0D0*DVZ2*N/V

         TORSHELL(3*D-1,3*B)= TORSHELL(3*D-1,3*B)+ M*DNZ2 + N*DMZ2
         TORSHELL(3*B,3*D-1)= TORSHELL(3*B,3*D-1)+ M*DNZ2 + N*DMZ2
C
C DY4DX2,DX2DY4
C
         N=(V*DUY4 - U*DVY4)/V**2
         VIXEN4= -2.0D0*(Y(B)-Y(C))*(X(B)-X(C))/RBC4
         VIXEN5= YPP*YE*VIXEN4 + (1.0D0+YE*DMUY4)*YE*DMUX2 + XPP*(XE*VIXEN4-DMUY4)
     1          + XE*DMUY4*(XE*DMUX2-1.0D0-MU) + ZPP*ZE*VIXEN4 + ZE**2*DMUY4*DMUX2 - DQY4*DQX2
         VIXEN3=(VIXEN1/VIXEN2)*VIXEN5 + DPX2*DQY4
         VIXEN6= YP*YE*VIXEN4 + (1.0D0+YE*DMUY4)*YE*DLAMBDAX2 + XP*(XE*VIXEN4-DMUY4)
     1          + XE*DMUY4*(XE*DLAMBDAX2-1.0D0-LAMBDA) + ZP*ZE*VIXEN4 + ZE**2*DMUY4*DLAMBDAX2
         DNX2= (V*VIXEN6 + DVX2*DUY4 - U*VIXEN3 - DUX2*DVY4)/V**2 - 2.0D0*DVX2*N/V

         TORSHELL(3*D-1,3*B-2)= TORSHELL(3*D-1,3*B-2)+ M*DNX2 + N*DMX2
         TORSHELL(3*B-2,3*D-1)= TORSHELL(3*B-2,3*D-1)+ M*DNX2 + N*DMX2
C
C DZ4DX2,DX2DZ4
C
         N=(V*DUZ4 - U*DVZ4)/V**2
         VIXEN4= -2.0D0*(Z(B)-Z(C))*(X(B)-X(C))/RBC4
         VIXEN5= ZPP*ZE*VIXEN4 + (1.0D0+ZE*DMUZ4)*ZE*DMUX2 + XPP*(XE*VIXEN4-DMUZ4)
     1          + XE*DMUZ4*(XE*DMUX2-1.0D0-MU) + YPP*YE*VIXEN4 + YE**2*DMUZ4*DMUX2 - DQZ4*DQX2
         VIXEN3=(VIXEN1/VIXEN2)*VIXEN5 + DPX2*DQZ4
         VIXEN6= ZP*ZE*VIXEN4 + (1.0D0+ZE*DMUZ4)*ZE*DLAMBDAX2 + XP*(XE*VIXEN4-DMUZ4)
     1          + XE*DMUZ4*(XE*DLAMBDAX2-1.0D0-LAMBDA) + YP*YE*VIXEN4 + YE**2*DMUZ4*DLAMBDAX2
         DNX2= (V*VIXEN6 + DVX2*DUZ4 - U*VIXEN3 - DUX2*DVZ4)/V**2 - 2.0D0*DVX2*N/V

         TORSHELL(3*D,3*B-2)= TORSHELL(3*D,3*B-2)+ M*DNX2 + N*DMX2
         TORSHELL(3*B-2,3*D)= TORSHELL(3*B-2,3*D)+ M*DNX2 + N*DMX2
C
C DZ4DY2,DY2DZ4
C
         N=(V*DUZ4 - U*DVZ4)/V**2
         VIXEN4= -2.0D0*(Z(B)-Z(C))*(Y(B)-Y(C))/RBC4
         VIXEN5= ZPP*ZE*VIXEN4 + (1.0D0+ZE*DMUZ4)*ZE*DMUY2 + YPP*(YE*VIXEN4-DMUZ4)
     1          + YE*DMUZ4*(YE*DMUY2-1.0D0-MU) + XPP*XE*VIXEN4 + XE**2*DMUZ4*DMUY2 - DQZ4*DQY2
         VIXEN3=(VIXEN1/VIXEN2)*VIXEN5 + DPY2*DQZ4
         VIXEN6= ZP*ZE*VIXEN4 + (1.0D0+ZE*DMUZ4)*ZE*DLAMBDAY2 + YP*(YE*VIXEN4-DMUZ4)
     1          + YE*DMUZ4*(YE*DLAMBDAY2-1.0D0-LAMBDA) + XP*XE*VIXEN4 + XE**2*DMUZ4*DLAMBDAY2
         DNY2= (V*VIXEN6 + DVY2*DUZ4 - U*VIXEN3 - DUY2*DVZ4)/V**2 - 2.0D0*DVY2*N/V

         TORSHELL(3*D,3*B-1)= TORSHELL(3*D,3*B-1)+ M*DNY2 + N*DMY2
         TORSHELL(3*B-1,3*D)= TORSHELL(3*B-1,3*D)+ M*DNY2 + N*DMY2
C
C DX4DX3,DX3DX4
C
         N=(V*DUX4 - U*DVX4)/V**2
         VIXEN4= 2.0D0*(X(B)-X(C))**2/RBC4 - 1.0D0/RBC2
         VIXEN5= XPP*(XE*VIXEN4+DMUX4) + (1.0D0+XE*DMUX4)*(XE*DMUX3+MU) + YPP*YE*VIXEN4
     1          + YE**2*DMUX4*DMUX3 + ZPP*ZE*VIXEN4 + ZE**2*DMUX4*DMUX3 - DQX4*DQX3
         VIXEN3=(VIXEN1/VIXEN2)*VIXEN5 + DPX3*DQX4
         VIXEN6= XP*(XE*VIXEN4+DMUX4) + (1.0D0+XE*DMUX4)*(XE*DLAMBDAX3+LAMBDA) + VIXEN4*
     1          (YP*YE+ZP*ZE) + (YE**2+ZE**2)*DMUX4*DLAMBDAX3
         DNX3= (V*VIXEN6 + DVX3*DUX4 - U*VIXEN3 - DUX3*DVX4)/V**2 - 2.0D0*DVX3*N/V

         TORSHELL(3*D-2,3*C-2)= TORSHELL(3*D-2,3*C-2)+ M*DNX3 + N*DMX3
         TORSHELL(3*C-2,3*D-2)= TORSHELL(3*C-2,3*D-2)+ M*DNX3 + N*DMX3
C
C DY4DY3,DY3DY4
C
         N=(V*DUY4 - U*DVY4)/V**2
         VIXEN4= 2.0D0*(Y(B)-Y(C))**2/RBC4 - 1.0D0/RBC2
         VIXEN5= YPP*(YE*VIXEN4+DMUY4) + (1.0D0+YE*DMUY4)*(YE*DMUY3+MU) + XPP*XE*VIXEN4
     1          + XE**2*DMUY4*DMUY3 + ZPP*ZE*VIXEN4 + ZE**2*DMUY4*DMUY3 - DQY4*DQY3
         VIXEN3=(VIXEN1/VIXEN2)*VIXEN5 + DPY3*DQY4
         VIXEN6= YP*(YE*VIXEN4+DMUY4) + (1.0D0+YE*DMUY4)*(YE*DLAMBDAY3+LAMBDA) + VIXEN4*
     1          (XP*XE+ZP*ZE) + (XE**2+ZE**2)*DMUY4*DLAMBDAY3
         DNY3= (V*VIXEN6 + DVY3*DUY4 - U*VIXEN3 - DUY3*DVY4)/V**2 - 2.0D0*DVY3*N/V

         TORSHELL(3*D-1,3*C-1)= TORSHELL(3*D-1,3*C-1)+ M*DNY3 + N*DMY3
         TORSHELL(3*C-1,3*D-1)= TORSHELL(3*C-1,3*D-1)+ M*DNY3 + N*DMY3
C
C DZ4DZ3,DZ3DZ4
C
         N=(V*DUZ4 - U*DVZ4)/V**2
         VIXEN4= 2.0D0*(Z(B)-Z(C))**2/RBC4 - 1.0D0/RBC2
         VIXEN5= ZPP*(ZE*VIXEN4+DMUZ4) + (1.0D0+ZE*DMUZ4)*(ZE*DMUZ3+MU) + XPP*XE*VIXEN4
     1          + XE**2*DMUZ4*DMUZ3 + YPP*YE*VIXEN4 + YE**2*DMUZ4*DMUZ3 - DQZ4*DQZ3
         VIXEN3=(VIXEN1/VIXEN2)*VIXEN5 + DPZ3*DQZ4
         VIXEN6= ZP*(ZE*VIXEN4+DMUZ4) + (1.0D0+ZE*DMUZ4)*(ZE*DLAMBDAZ3+LAMBDA) + VIXEN4*
     1          (XP*XE+YP*YE) + (XE**2+YE**2)*DMUZ4*DLAMBDAZ3
         DNZ3= (V*VIXEN6 + DVZ3*DUZ4 - U*VIXEN3 - DUZ3*DVZ4)/V**2 - 2.0D0*DVZ3*N/V

         TORSHELL(3*D,3*C)= TORSHELL(3*D,3*C)+ M*DNZ3 + N*DMZ3
         TORSHELL(3*C,3*D)= TORSHELL(3*C,3*D)+ M*DNZ3 + N*DMZ3
C
C DX4DY3,DY3DX4
C
         N=(V*DUX4 - U*DVX4)/V**2
         VIXEN4= 2.0D0*(X(B)-X(C))*(Y(B)-Y(C))/RBC4
         VIXEN5= XPP*XE*VIXEN4 + (1.0D0+XE*DMUX4)*XE*DMUY3 + YPP*(YE*VIXEN4+DMUX4)
     1          + YE*DMUX4*(YE*DMUY3+MU) + ZPP*ZE*VIXEN4 + ZE**2*DMUX4*DMUY3 - DQX4*DQY3
         VIXEN3=(VIXEN1/VIXEN2)*VIXEN5 + DPY3*DQX4
         VIXEN6= XP*XE*VIXEN4 + (1.0D0+XE*DMUX4)*XE*DLAMBDAY3 + YP*(YE*VIXEN4+DMUX4)
     1          + YE*DMUX4*(YE*DLAMBDAY3+LAMBDA) + ZP*ZE*VIXEN4 + ZE**2*DMUX4*DLAMBDAY3
         DNY3= (V*VIXEN6 + DVY3*DUX4 - U*VIXEN3 - DUY3*DVX4)/V**2 - 2.0D0*DVY3*N/V

         TORSHELL(3*D-2,3*C-1)= TORSHELL(3*D-2,3*C-1) + M*DNY3 + N*DMY3
         TORSHELL(3*C-1,3*D-2)= TORSHELL(3*C-1,3*D-2) + M*DNY3 + N*DMY3
C
C DX4DZ3,DZ3DX4
C
         N=(V*DUX4 - U*DVX4)/V**2
         VIXEN4= 2.0D0*(X(B)-X(C))*(Z(B)-Z(C))/RBC4
         VIXEN5= XPP*XE*VIXEN4 + (1.0D0+XE*DMUX4)*XE*DMUZ3 + ZPP*(ZE*VIXEN4+DMUX4)
     1          + ZE*DMUX4*(ZE*DMUZ3+MU) + YPP*YE*VIXEN4 + YE**2*DMUX4*DMUZ3 - DQX4*DQZ3
         VIXEN3=(VIXEN1/VIXEN2)*VIXEN5 + DPZ3*DQX4
         VIXEN6= XP*XE*VIXEN4 + (1.0D0+XE*DMUX4)*XE*DLAMBDAZ3 + ZP*(ZE*VIXEN4+DMUX4)
     1          + ZE*DMUX4*(ZE*DLAMBDAZ3+LAMBDA) + YP*YE*VIXEN4 + YE**2*DMUX4*DLAMBDAZ3
         DNZ3= (V*VIXEN6 + DVZ3*DUX4 - U*VIXEN3 - DUZ3*DVX4)/V**2 - 2.0D0*DVZ3*N/V

         TORSHELL(3*D-2,3*C)= TORSHELL(3*D-2,3*C) + M*DNZ3 + N*DMZ3
         TORSHELL(3*C,3*D-2)= TORSHELL(3*C,3*D-2) + M*DNZ3 + N*DMZ3
C
C DY4DZ3,DZ3DY4
C
         N=(V*DUY4 - U*DVY4)/V**2
         VIXEN4= 2.0D0*(Y(B)-Y(C))*(Z(B)-Z(C))/RBC4
         VIXEN5= YPP*YE*VIXEN4 + (1.0D0+YE*DMUY4)*YE*DMUZ3 + ZPP*(ZE*VIXEN4+DMUY4)
     1          + ZE*DMUY4*(ZE*DMUZ3+MU) + XPP*XE*VIXEN4 + XE**2*DMUY4*DMUZ3 - DQY4*DQZ3
         VIXEN3=(VIXEN1/VIXEN2)*VIXEN5 + DPZ3*DQY4
         VIXEN6= YP*YE*VIXEN4 + (1.0D0+YE*DMUY4)*YE*DLAMBDAZ3 + ZP*(ZE*VIXEN4+DMUY4)
     1          + ZE*DMUY4*(ZE*DLAMBDAZ3+LAMBDA) + XP*XE*VIXEN4 + XE**2*DMUY4*DLAMBDAZ3
         DNZ3= (V*VIXEN6 + DVZ3*DUY4 - U*VIXEN3 - DUZ3*DVY4)/V**2 - 2.0D0*DVZ3*N/V

         TORSHELL(3*D-1,3*C)= TORSHELL(3*D-1,3*C) + M*DNZ3 + N*DMZ3
         TORSHELL(3*C,3*D-1)= TORSHELL(3*C,3*D-1) + M*DNZ3 + N*DMZ3
C
C DY4DX3,DX3DY4
C
         N=(V*DUY4 - U*DVY4)/V**2
         VIXEN4= 2.0D0*(Y(B)-Y(C))*(X(B)-X(C))/RBC4
         VIXEN5= YPP*YE*VIXEN4 + (1.0D0+YE*DMUY4)*YE*DMUX3 + XPP*(XE*VIXEN4+DMUY4)
     1          + XE*DMUY4*(XE*DMUX3+MU) + ZPP*ZE*VIXEN4 + ZE**2*DMUY4*DMUX3 - DQY4*DQX3
         VIXEN3=(VIXEN1/VIXEN2)*VIXEN5 + DPX3*DQY4
         VIXEN6= YP*YE*VIXEN4 + (1.0D0+YE*DMUY4)*YE*DLAMBDAX3 + XP*(XE*VIXEN4+DMUY4)
     1          + XE*DMUY4*(XE*DLAMBDAX3+LAMBDA) + ZP*ZE*VIXEN4 + ZE**2*DMUY4*DLAMBDAX3
         DNX3= (V*VIXEN6 + DVX3*DUY4 - U*VIXEN3 - DUX3*DVY4)/V**2 - 2.0D0*DVX3*N/V

         TORSHELL(3*D-1,3*C-2)= TORSHELL(3*D-1,3*C-2) + M*DNX3 + N*DMX3
         TORSHELL(3*C-2,3*D-1)= TORSHELL(3*C-2,3*D-1) + M*DNX3 + N*DMX3
C
C DZ4DX3,DX3DZ4
C
         N=(V*DUZ4 - U*DVZ4)/V**2
         VIXEN4= 2.0D0*(Z(B)-Z(C))*(X(B)-X(C))/RBC4
         VIXEN5= ZPP*ZE*VIXEN4 + (1.0D0+ZE*DMUZ4)*ZE*DMUX3 + XPP*(XE*VIXEN4+DMUZ4)
     1          + XE*DMUZ4*(XE*DMUX3+MU) + YPP*YE*VIXEN4 + YE**2*DMUZ4*DMUX3 - DQZ4*DQX3
         VIXEN3=(VIXEN1/VIXEN2)*VIXEN5 + DPX3*DQZ4
         VIXEN6= ZP*ZE*VIXEN4 + (1.0D0+ZE*DMUZ4)*ZE*DLAMBDAX3 + XP*(XE*VIXEN4+DMUZ4)
     1          + XE*DMUZ4*(XE*DLAMBDAX3+LAMBDA) + YP*YE*VIXEN4 + YE**2*DMUZ4*DLAMBDAX3
         DNX3= (V*VIXEN6 + DVX3*DUZ4 - U*VIXEN3 - DUX3*DVZ4)/V**2 - 2.0D0*DVX3*N/V

         TORSHELL(3*D,3*C-2)= TORSHELL(3*D,3*C-2) + M*DNX3 + N*DMX3
         TORSHELL(3*C-2,3*D)= TORSHELL(3*C-2,3*D) + M*DNX3 + N*DMX3
C
C DZ4DY3,DY3DZ4
C
         N=(V*DUZ4 - U*DVZ4)/V**2
         VIXEN4= 2.0D0*(Z(B)-Z(C))*(Y(B)-Y(C))/RBC4
         VIXEN5= ZPP*ZE*VIXEN4 + (1.0D0+ZE*DMUZ4)*ZE*DMUY3 + YPP*(YE*VIXEN4+DMUZ4)
     1          + YE*DMUZ4*(YE*DMUY3+MU) + XPP*XE*VIXEN4 + XE**2*DMUZ4*DMUY3 - DQZ4*DQY3
         VIXEN3=(VIXEN1/VIXEN2)*VIXEN5 + DPY3*DQZ4
         VIXEN6= ZP*ZE*VIXEN4 + (1.0D0+ZE*DMUZ4)*ZE*DLAMBDAY3 + YP*(YE*VIXEN4+DMUZ4)
     1          + YE*DMUZ4*(YE*DLAMBDAY3+LAMBDA) + XP*XE*VIXEN4 + XE**2*DMUZ4*DLAMBDAY3
         DNY3= (V*VIXEN6 + DVY3*DUZ4 - U*VIXEN3 - DUY3*DVZ4)/V**2 - 2.0D0*DVY3*N/V

         TORSHELL(3*D,3*C-1)= TORSHELL(3*D,3*C-1) + M*DNY3 + N*DMY3
         TORSHELL(3*C-1,3*D)= TORSHELL(3*C-1,3*D) + M*DNY3 + N*DMY3
C
C DX3DX3
C
         N=(V*DUX3 - U*DVX3)/V**2
         VIXEN4= 8.0D0*(XE**2-RBC2)*MUCRAP/RBC6 - 4.0D0*XE*(X(B)-X(D))/RBC4
         VIXEN5= 8.0D0*(XE**2-RBC2)*LAMBDACRAP/RBC6 - 4.0D0*XE*(X(B)-X(A))/RBC4
         VIXEN7= XPP*(2.0D0*DMUX3+XE*VIXEN4) + (MU+XE*DMUX3)**2 + VIXEN4*(YPP*YE+ZPP*ZE)
     1          + DMUX3**2*(YE**2+ZE**2) - DQX3**2
         VIXEN8= XP*(2.0D0*DLAMBDAX3+XE*VIXEN5) + (LAMBDA+XE*DLAMBDAX3)**2 + VIXEN5*(YP*YE+ZP*ZE)
     1          + DLAMBDAX3**2*(YE**2+ZE**2) - DPX3**2
         VIXEN3= (VIXEN1/VIXEN2)*VIXEN7 + (VIXEN2/VIXEN1)*VIXEN8 + 2.0D0*DPX3*DQX3
         VIXEN6= XP*(2.0D0*DMUX3+XE*VIXEN4) + 2.0D0*(MU+XE*DMUX3)*(LAMBDA+XE*DLAMBDAX3) 
     1          + XPP*(2.0D0*DLAMBDAX3+XE*VIXEN5) + VIXEN4*(YP*YE+ZP*ZE) 
     2          + VIXEN5*(YPP*YE+ZPP*ZE) + 2.0D0*DLAMBDAX3*DMUX3*(YE**2+ZE**2)
         DNX3= (V*VIXEN6 - U*VIXEN3)/V**2 - 2.0D0*DVX3*N/V

         TORSHELL(3*C-2,3*C-2)= TORSHELL(3*C-2,3*C-2)+ M*DNX3 + N*DMX3
C
C DY3DY3
C
         N=(V*DUY3 - U*DVY3)/V**2
         VIXEN4= 8.0D0*(XE**2-RBC2)*MUCRAP/RBC6 - 4.0D0*XE*(X(B)-X(D))/RBC4
         VIXEN5= 8.0D0*(XE**2-RBC2)*LAMBDACRAP/RBC6 - 4.0D0*XE*(X(B)-X(A))/RBC4
         VIXEN7= YPP*(2.0D0*DMUY3+YE*VIXEN4) + (MU+YE*DMUY3)**2 + VIXEN4*(XPP*XE+ZPP*ZE)
     1          + DMUY3**2*(XE**2+ZE**2) - DQY3**2
         VIXEN8= YP*(2.0D0*DLAMBDAY3+YE*VIXEN5) + (LAMBDA+YE*DLAMBDAY3)**2 + VIXEN5*(XP*XE+ZP*ZE)
     1          + DLAMBDAY3**2*(XE**2+ZE**2) - DPY3**2
         VIXEN3= (VIXEN1/VIXEN2)*VIXEN7 + (VIXEN2/VIXEN1)*VIXEN8 + 2.0D0*DPY3*DQY3
         VIXEN6= YP*(2.0D0*DMUY3+YE*VIXEN4) + 2.0D0*(MU+YE*DMUY3)*(LAMBDA+YE*DLAMBDAY3)
     1          + YPP*(2.0D0*DLAMBDAY3+YE*VIXEN5) + VIXEN4*(XP*XE+ZP*ZE)
     2          + VIXEN5*(XPP*XE+ZPP*ZE) + 2.0D0*DLAMBDAY3*DMUY3*(XE**2+ZE**2)
         DNY3= (V*VIXEN6 - U*VIXEN3)/V**2 - 2.0D0*DVY3*N/V

         TORSHELL(3*C-1,3*C-1)= TORSHELL(3*C-1,3*C-1)+ M*DNY3 + N*DMY3
C
C DZ3DZ3
C
         N=(V*DUZ3 - U*DVZ3)/V**2
         VIXEN4= 8.0D0*(XE**2-RBC2)*MUCRAP/RBC6 - 4.0D0*XE*(X(B)-X(D))/RBC4
         VIXEN5= 8.0D0*(XE**2-RBC2)*LAMBDACRAP/RBC6 - 4.0D0*XE*(X(B)-X(A))/RBC4
         VIXEN7= ZPP*(2.0D0*DMUZ3+ZE*VIXEN4) + (MU+ZE*DMUZ3)**2 + VIXEN4*(XPP*XE+YPP*YE)
     1          + DMUZ3**2*(XE**2+YE**2) - DQZ3**2
         VIXEN8= ZP*(2.0D0*DLAMBDAZ3+ZE*VIXEN5) + (LAMBDA+ZE*DLAMBDAZ3)**2 + VIXEN5*(XP*XE+YP*YE)
     1          + DLAMBDAZ3**2*(XE**2+YE**2) - DPZ3**2
         VIXEN3= (VIXEN1/VIXEN2)*VIXEN7 + (VIXEN2/VIXEN1)*VIXEN8 + 2.0D0*DPZ3*DQZ3
         VIXEN6= ZP*(2.0D0*DMUZ3+ZE*VIXEN4) + 2.0D0*(MU+ZE*DMUZ3)*(LAMBDA+ZE*DLAMBDAZ3)
     1          + ZPP*(2.0D0*DLAMBDAZ3+ZE*VIXEN5) + VIXEN4*(XP*XE+YP*YE)
     2          + VIXEN5*(XPP*XE+YPP*YE) + 2.0D0*DLAMBDAZ3*DMUZ3*(XE**2+YE**2)
         DNZ3= (V*VIXEN6 - U*VIXEN3)/V**2 - 2.0D0*DVZ3*N/V

         TORSHELL(3*C,3*C)= TORSHELL(3*C,3*C)+ M*DNZ3 + N*DMZ3
C
C DX3DY3,DY3DX3
C
         N=(V*DUX3 - U*DVX3)/V**2
         VIXEN4= 8.0D0*MUCRAP*XE*YE/RBC6 - 2.0D0*(YE*(X(B)-X(D))+XE*(Y(B)-Y(D)))/RBC4
         VIXEN5= 8.0D0*MUCRAP*XE*YE/RBC6 - 2.0D0*(YE*(X(B)-X(A))+XE*(Y(B)-Y(A)))/RBC4
         VIXEN7= XPP*(DMUY3+XE*VIXEN4) + (MU+XE*DMUX3)*XE*DMUY3 + YPP*(YE*VIXEN4+DMUX3)
     1          + YE*DMUX3*(MU+YE*DMUY3) + ZPP*ZE*VIXEN4 + ZE**2*DMUX3*DMUY3 - DQX3*DQY3
         VIXEN8= XP*(DLAMBDAY3+XE*VIXEN5) + (LAMBDA+XE*DLAMBDAX3)*XE*DLAMBDAY3 
     1          + YP*(YE*VIXEN5 + DLAMBDAX3) + YE*DLAMBDAX3*(LAMBDA+YE*DLAMBDAY3) 
     2          + ZP*ZE*VIXEN5 + ZE**2*DLAMBDAX3*DLAMBDAY3 - DPX3*DPY3
         VIXEN3= (VIXEN1/VIXEN2)*VIXEN7 + (VIXEN2/VIXEN1)*VIXEN8 + DPX3*DQY3 + DPY3*DQX3
         VIXEN6= XP*(DMUY3+XE*VIXEN4) + (MU+XE*DMUX3)*XE*DLAMBDAY3 + XPP*(DLAMBDAY3+XE*VIXEN5)
     1          + (LAMBDA+XE*DLAMBDAX3)*XE*DMUY3 + YP*(YE*VIXEN4+DMUX3) 
     2          + YE*DMUX3*(YE*DLAMBDAY3+LAMBDA) + YPP*(YE*VIXEN5 + DLAMBDAX3) 
     3          + YE*DLAMBDAX3*(MU+YE*DMUY3) + ZP*ZE*VIXEN4 + ZPP*ZE*VIXEN5
     4          + ZE**2*(DLAMBDAX3*DMUY3+DMUX3*DLAMBDAY3)
         DNY3= (V*VIXEN6 + DVY3*DUX3 - U*VIXEN3 - DUY3*DVX3)/V**2 - 2.0D0*DVY3*N/V

         TORSHELL(3*C-2,3*C-1)= TORSHELL(3*C-2,3*C-1)+ M*DNY3 + N*DMY3
         TORSHELL(3*C-1,3*C-2)= TORSHELL(3*C-1,3*C-2)+ M*DNY3 + N*DMY3
C
C DX3DZ3,DZ3DX3
C
         N=(V*DUX3 - U*DVX3)/V**2
         VIXEN4= 8.0D0*MUCRAP*XE*ZE/RBC6 - 2.0D0*(ZE*(X(B)-X(D))+XE*(Z(B)-Z(D)))/RBC4
         VIXEN5= 8.0D0*MUCRAP*XE*ZE/RBC6 - 2.0D0*(ZE*(X(B)-X(A))+XE*(Z(B)-Z(A)))/RBC4
         VIXEN7= XPP*(DMUZ3+XE*VIXEN4) + (MU+XE*DMUX3)*XE*DMUZ3 + ZPP*(ZE*VIXEN4+DMUX3)
     1          + ZE*DMUX3*(MU+ZE*DMUZ3) + YPP*YE*VIXEN4 + YE**2*DMUX3*DMUZ3 - DQX3*DQZ3
         VIXEN8= XP*(DLAMBDAZ3+XE*VIXEN5) + (LAMBDA+XE*DLAMBDAX3)*XE*DLAMBDAZ3
     1          + ZP*(ZE*VIXEN5 + DLAMBDAX3) + ZE*DLAMBDAX3*(LAMBDA+ZE*DLAMBDAZ3)
     2          + YP*YE*VIXEN5 + YE**2*DLAMBDAX3*DLAMBDAZ3 - DPX3*DPZ3
         VIXEN3= (VIXEN1/VIXEN2)*VIXEN7 + (VIXEN2/VIXEN1)*VIXEN8 + DPX3*DQZ3 + DPZ3*DQX3
         VIXEN6= XP*(DMUZ3+XE*VIXEN4) + (MU+XE*DMUX3)*XE*DLAMBDAZ3 + XPP*(DLAMBDAZ3+XE*VIXEN5)
     1          + (LAMBDA+XE*DLAMBDAX3)*XE*DMUZ3 + ZP*(ZE*VIXEN4+DMUX3)
     2          + ZE*DMUX3*(ZE*DLAMBDAZ3+LAMBDA) + ZPP*(ZE*VIXEN5 + DLAMBDAX3)
     3          + ZE*DLAMBDAX3*(MU+ZE*DMUZ3) + YP*YE*VIXEN4 + YPP*YE*VIXEN5
     4          + YE**2*(DLAMBDAX3*DMUZ3+DMUX3*DLAMBDAZ3)
         DNZ3= (V*VIXEN6 + DVZ3*DUX3 - U*VIXEN3 - DUZ3*DVX3)/V**2 - 2.0D0*DVZ3*N/V

         TORSHELL(3*C-2,3*C)= TORSHELL(3*C-2,3*C)+ M*DNZ3 + N*DMZ3
         TORSHELL(3*C,3*C-2)= TORSHELL(3*C,3*C-2)+ M*DNZ3 + N*DMZ3
C
C DY3DZ3,DZ3DY3
C
         N=(V*DUY3 - U*DVY3)/V**2
         VIXEN4= 8.0D0*MUCRAP*YE*ZE/RBC6 - 2.0D0*(ZE*(Y(B)-Y(D))+YE*(Z(B)-Z(D)))/RBC4
         VIXEN5= 8.0D0*MUCRAP*YE*ZE/RBC6 - 2.0D0*(ZE*(Y(B)-Y(A))+YE*(Z(B)-Z(A)))/RBC4
         VIXEN7= YPP*(DMUZ3+YE*VIXEN4) + (MU+YE*DMUY3)*YE*DMUZ3 + ZPP*(ZE*VIXEN4+DMUY3)
     1          + ZE*DMUY3*(MU+ZE*DMUZ3) + XPP*XE*VIXEN4 + XE**2*DMUY3*DMUZ3 - DQY3*DQZ3
         VIXEN8= YP*(DLAMBDAZ3+YE*VIXEN5) + (LAMBDA+YE*DLAMBDAY3)*YE*DLAMBDAZ3
     1          + ZP*(ZE*VIXEN5 + DLAMBDAY3) + ZE*DLAMBDAY3*(LAMBDA+ZE*DLAMBDAZ3)
     2          + XP*XE*VIXEN5 + XE**2*DLAMBDAY3*DLAMBDAZ3 - DPY3*DPZ3
         VIXEN3= (VIXEN1/VIXEN2)*VIXEN7 + (VIXEN2/VIXEN1)*VIXEN8 + DPY3*DQZ3 + DPZ3*DQY3
         VIXEN6= YP*(DMUZ3+YE*VIXEN4) + (MU+YE*DMUY3)*YE*DLAMBDAZ3 + YPP*(DLAMBDAZ3+YE*VIXEN5)
     1          + (LAMBDA+YE*DLAMBDAY3)*YE*DMUZ3 + ZP*(ZE*VIXEN4+DMUY3)
     2          + ZE*DMUY3*(ZE*DLAMBDAZ3+LAMBDA) + ZPP*(ZE*VIXEN5 + DLAMBDAY3)
     3          + ZE*DLAMBDAY3*(MU+ZE*DMUZ3) + XP*XE*VIXEN4 + XPP*XE*VIXEN5
     4          + XE**2*(DLAMBDAY3*DMUZ3+DMUY3*DLAMBDAZ3)
         DNZ3= (V*VIXEN6 + DVZ3*DUY3 - U*VIXEN3 - DUZ3*DVY3)/V**2 - 2.0D0*DVZ3*N/V

         TORSHELL(3*C-1,3*C)= TORSHELL(3*C-1,3*C)+ M*DNZ3 + N*DMZ3
         TORSHELL(3*C,3*C-1)= TORSHELL(3*C,3*C-1)+ M*DNZ3 + N*DMZ3
C
C DX2DX2
C
         N=(V*DUX2 - U*DVX2)/V**2
         VIXEN4= -2.0D0/RBC2 - (4.0D0*(X(B)-X(C))*(X(C)+X(D)+2.0D0*X(B))-
     1          2.0D0*MUCRAP)/RBC4 + 8.0D0*(X(B)-X(C))**2*MUCRAP/RBC6
         VIXEN5= -2.0D0/RBC2 - (4.0D0*(X(B)-X(C))*(X(C)+X(A)+2.0D0*X(B))-
     1    2.0D0*LAMBDACRAP)/RBC4 + 8.0D0*(X(B)-X(C))**2*LAMBDACRAP/RBC6
         VIXEN7= XP*(XE*VIXEN5-2.0D0*DLAMBDAX2) + (XE*DLAMBDAX2
     1    -1.0D0-LAMBDA)**2 + (YP*YE+ZP*ZE)*VIXEN5 - DPX2**2
     2    + (YE**2+ZE**2)*DLAMBDAX2**2
         VIXEN8= XPP*(XE*VIXEN4-2.0D0*DMUX2) + (XE*DMUX2-1.0D0-MU)
     1    **2 + (YPP*YE+ZPP*ZE)*VIXEN4 - DQX2**2
     2    + (YE**2+ZE**2)*DMUX2**2
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN7 + (VIXEN1/VIXEN2)*VIXEN8+2.0D0*DPX2*DQX2
         VIXEN6= XP*(XE*VIXEN4-2.0D0*DMUX2) + 2.0D0*(XE*DMUX2-1-MU)
     1    *(XE*DLAMBDAX2-1.0D0-LAMBDA) + XPP*(XE*VIXEN5-2.0D0*
     2    DLAMBDAX2) + VIXEN4*(YP*YE+ZP*ZE) + VIXEN5*
     3    (YPP*YE+ZPP*ZE) + 2.0D0*DMUX2*DLAMBDAX2*(YE
     4    **2+ZE**2)
         DNX2= (V*VIXEN6 - U*VIXEN3)/V**2 - 2.0D0*DVX2*N/V

         TORSHELL(3*B-2,3*B-2)= TORSHELL(3*B-2,3*B-2)+ M*DNX2 + N*DMX2
C
C DY2DY2
C
         N=(V*DUY2 - U*DVY2)/V**2
         VIXEN4= -2.0D0/RBC2 - (4.0D0*(Y(B)-Y(C))*(Y(C)+Y(D)+2.0D0*Y(B))-
     1          2.0D0*MUCRAP)/RBC4 + 8.0D0*(Y(B)-Y(C))**2*MUCRAP/RBC6
         VIXEN5= -2.0D0/RBC2 - (4.0D0*(Y(B)-Y(C))*(Y(C)+Y(A)+2.0D0*Y(B))-
     1    2.0D0*LAMBDACRAP)/RBC4 + 8.0D0*(Y(B)-Y(C))**2*LAMBDACRAP/RBC6
         VIXEN7= YP*(YE*VIXEN5-2.0D0*DLAMBDAY2) + (YE*DLAMBDAY2
     1    -1.0D0-LAMBDA)**2 + (XP*XE+ZP*ZE)*VIXEN5 - DPY2**2
     2    + (XE**2+ZE**2)*DLAMBDAY2**2
         VIXEN8= YPP*(YE*VIXEN4-2.0D0*DMUY2) + (YE*DMUY2-1.0D0-MU)
     1    **2 + (XPP*XE+ZPP*ZE)*VIXEN4 - DQY2**2
     2    + (XE**2+ZE**2)*DMUY2**2
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN7 + (VIXEN1/VIXEN2)*VIXEN8+2.0D0*DPY2*DQY2
         VIXEN6= YP*(YE*VIXEN4-2.0D0*DMUY2) + 2.0D0*(YE*DMUY2-1.0D0-MU)
     1    *(YE*DLAMBDAY2-1.0D0-LAMBDA) + YPP*(YE*VIXEN5-2.0D0*
     2    DLAMBDAY2) + VIXEN4*(XP*XE+ZP*ZE) + VIXEN5*
     3    (XPP*XE+ZPP*ZE) + 2.0D0*DMUY2*DLAMBDAY2*(XE
     4    **2+ZE**2)
         DNY2= (V*VIXEN6 - U*VIXEN3)/V**2 - 2.0D0*DVY2*N/V

         TORSHELL(3*B-1,3*B-1)= TORSHELL(3*B-1,3*B-1)+ M*DNY2 + N*DMY2
C
C DZ2DZ2
C
         N=(V*DUZ2 - U*DVZ2)/V**2
         VIXEN4= -2.0D0/RBC2 - (4.0D0*(Z(B)-Z(C))*(Z(C)+Z(D)+2.0D0*Z(B))-
     1          2.0D0*MUCRAP)/RBC4 + 8.0D0*(Z(B)-Z(C))**2*MUCRAP/RBC6
         VIXEN5= -2.0D0/RBC2 - (4.0D0*(Z(B)-Z(C))*(Z(C)+Z(A)+2.0D0*Z(B))-
     1    2.0D0*LAMBDACRAP)/RBC4 + 8.0D0*(Z(B)-Z(C))**2*LAMBDACRAP/RBC6
         VIXEN7= ZP*(ZE*VIXEN5-2.0D0*DLAMBDAZ2) + (ZE*DLAMBDAZ2
     1    -1.0D0-LAMBDA)**2 + (XP*XE+YP*YE)*VIXEN5 - DPZ2**2
     2    + (XE**2+YE**2)*DLAMBDAZ2**2
         VIXEN8= ZPP*(ZE*VIXEN4-2.0D0*DMUZ2) + (ZE*DMUZ2-1.0D0-MU)
     1    **2 + (XPP*XE+YPP*YE)*VIXEN4 - DQZ2**2
     2    + (XE**2+YE**2)*DMUZ2**2
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN7 + (VIXEN1/VIXEN2)*VIXEN8+2.0D0*DPZ2*DQZ2
         VIXEN6= ZP*(ZE*VIXEN4-2.0D0*DMUZ2) + 2.0D0*(ZE*DMUZ2-1.0D0-MU)
     1    *(ZE*DLAMBDAZ2-1.0D0-LAMBDA) + ZPP*(ZE*VIXEN5-2.0D0*
     2    DLAMBDAZ2) + VIXEN4*(XP*XE+YP*YE) + VIXEN5*
     3    (XPP*XE+YPP*YE) + 2.0D0*DMUZ2*DLAMBDAZ2*(XE
     4    **2+YE**2)
         DNZ2= (V*VIXEN6 - U*VIXEN3)/V**2 - 2.0D0*DVZ2*N/V

         TORSHELL(3*B,3*B)= TORSHELL(3*B,3*B)+ M*DNZ2 + N*DMZ2
C
C DX2DY2,DY2DX2
C
         N=(V*DUX2 - U*DVX2)/V**2
         VIXEN4= (-2.0D0*(Y(B)-Y(C))*(X(C)+X(D)-2.0D0*X(B))-2.0D0*(X(B)-X(C))*(Y(C)+Y(D)-2.0D0*Y(B)))/RBC4 
     1          + 8.0D0*(X(B)-X(C))*(Y(B)-Y(C))*MUCRAP
     2          /RBC6
         VIXEN5= (-2.0D0*(Y(B)-Y(C))*(X(C)+X(A)-2.0D0*X(B))-2.0D0*(X(B)-X(C))*(Y(C)+Y(A)-2.0D0*Y(B)))/RBC4 
     1          + 8.0D0*(X(B)-X(C))*(Y(B)-Y(C))
     2          *LAMBDACRAP/RBC6
         VIXEN7= XP*(XE*VIXEN5-DLAMBDAY2) + (XE*DLAMBDAX2-1.0D0-LAMBDA)*XE*DLAMBDAY2 + YP*(YE
     1          *VIXEN5-DLAMBDAX2) + YE*DLAMBDAX2*(YE*DLAMBDAY2-1.0D0-LAMBDA) + ZP*ZE*VIXEN5
     2          + ZE**2*DLAMBDAX2*DLAMBDAY2 - DPX2*DPY2
         VIXEN8= XPP*(XE*VIXEN4-DMUY2) + (XE*DMUX2-1.0D0-MU)*XE*DMUY2 + YPP*(YE
     1          *VIXEN4-DMUX2) + YE*DMUX2*(YE*DMUY2-1.0D0-MU) + ZPP*ZE*VIXEN4
     2          + ZE**2*DMUX2*DMUY2 - DQX2*DQY2
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN7 + (VIXEN1/VIXEN2)*VIXEN8 + DPX2*DQY2 + DQX2*DPY2
         VIXEN6= XP*(XE*VIXEN4-DMUY2) + (XE*DMUX2-1.0D0-MU)*XE*DLAMBDAY2 + XPP*(XE*VIXEN5-DLAMBDAY2)
     1          + (XE*DLAMBDAX2-1.0D0-LAMBDA)*XE*DMUY2 + YP*(YE*VIXEN4-DMUX2) + YE*DMUX2
     2          *(YE*DLAMBDAY2-1.0D0-LAMBDA) + YPP*(YE*VIXEN5-DLAMBDAX2) + YE*DLAMBDAX2
     3          *(YE*DMUY2-1.0D0-MU) + ZP*ZE*VIXEN4 + ZPP*ZE*VIXEN5 + ZE**2*(DMUX2*DLAMBDAY2
     4          +DMUY2*DLAMBDAX2)
         DNY2= (V*VIXEN6 + DVY2*DUX2 - U*VIXEN3 - DVX2*DUY2)/V**2 - 2.0D0*DVY2*N/V

         TORSHELL(3*B-2,3*B-1)= TORSHELL(3*B-2,3*B-1)+ M*DNY2 + N*DMY2
         TORSHELL(3*B-1,3*B-2)= TORSHELL(3*B-1,3*B-2)+ M*DNY2 + N*DMY2
C
C DX2DZ2,DZ2DX2
C
         N=(V*DUX2 - U*DVX2)/V**2
         VIXEN4= (-2.0D0*(Z(B)-Z(C))*(X(C)+X(D)-2.0D0*X(B))-2.0D0*(X(B)-X(C))*(Z(C)+Z(D)-2.0D0*Z(B)))/RBC4 
     1          + 8.0D0*(X(B)-X(C))*(Z(B)-Z(C))*MUCRAP
     2          /RBC6
         VIXEN5= (-2.0D0*(Z(B)-Z(C))*(X(C)+X(A)-2.0D0*X(B))-2.0D0*(X(B)-X(C))*(Z(C)+Z(A)-2.0D0*Z(B)))/RBC4 
     1          + 8.0D0*(X(B)-X(C))*(Z(B)-Z(C))
     2          *LAMBDACRAP/RBC6
         VIXEN7= XP*(XE*VIXEN5-DLAMBDAZ2) + (XE*DLAMBDAX2-1.0D0-LAMBDA)*XE*DLAMBDAZ2 + ZP*(ZE
     1          *VIXEN5-DLAMBDAX2) + ZE*DLAMBDAX2*(ZE*DLAMBDAZ2-1.0D0-LAMBDA) + YP*YE*VIXEN5
     2          + YE**2*DLAMBDAX2*DLAMBDAZ2 - DPX2*DPZ2
         VIXEN8= XPP*(XE*VIXEN4-DMUZ2) + (XE*DMUX2-1.0D0-MU)*XE*DMUZ2 + ZPP*(ZE
     1          *VIXEN4-DMUX2) + ZE*DMUX2*(ZE*DMUZ2-1.0D0-MU) + YPP*YE*VIXEN4
     2          + YE**2*DMUX2*DMUZ2 - DQX2*DQZ2
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN7 + (VIXEN1/VIXEN2)*VIXEN8 + DPX2*DQZ2 + DQX2*DPZ2
         VIXEN6= XP*(XE*VIXEN4-DMUZ2) + (XE*DMUX2-1.0D0-MU)*XE*DLAMBDAZ2 + XPP*(XE*VIXEN5-DLAMBDAZ2)
     1          + (XE*DLAMBDAX2-1.0D0-LAMBDA)*XE*DMUZ2 + ZP*(ZE*VIXEN4-DMUX2) + ZE*DMUX2
     2          *(ZE*DLAMBDAZ2-1.0D0-LAMBDA) + ZPP*(ZE*VIXEN5-DLAMBDAX2) + ZE*DLAMBDAX2
     3          *(ZE*DMUZ2-1.0D0-MU) + YP*YE*VIXEN4 + YPP*YE*VIXEN5 + YE**2*(DMUX2*DLAMBDAZ2
     4          +DMUZ2*DLAMBDAX2)
         DNZ2= (V*VIXEN6 + DVZ2*DUX2 - U*VIXEN3 - DVX2*DUZ2)/V**2 - 2.0D0*DVZ2*N/V

         TORSHELL(3*B-2,3*B)= TORSHELL(3*B-2,3*B)+ M*DNZ2 + N*DMZ2
         TORSHELL(3*B,3*B-2)= TORSHELL(3*B,3*B-2)+ M*DNZ2 + N*DMZ2
C
C DY2DZ2,DZ2DY2
C
         N=(V*DUY2 - U*DVY2)/V**2
         VIXEN4= (-2.0D0*(Z(B)-Z(C))*(Y(C)+Y(D)-2.0D0*Y(B))-2.0D0*(Y(B)-Y(C))*(Z(C)+Z(D)-2.0D0*Z(B)))/RBC4 
     1          + 8.0D0*(Y(B)-Y(C))*(Z(B)-Z(C))*MUCRAP
     2          /RBC6
         VIXEN5= (-2.0D0*(Z(B)-Z(C))*(Y(C)+Y(A)-2.0D0*Y(B))-2.0D0*(Y(B)-Y(C))*(Z(C)+Z(A)-2.0D0*Z(B)))/RBC4 
     1          + 8.0D0*(Y(B)-Y(C))*(Z(B)-Z(C))
     2          *LAMBDACRAP/RBC6
         VIXEN7= YP*(YE*VIXEN5-DLAMBDAZ2) + (YE*DLAMBDAY2-1.0D0-LAMBDA)*YE*DLAMBDAZ2 + ZP*(ZE
     1          *VIXEN5-DLAMBDAY2) + ZE*DLAMBDAY2*(ZE*DLAMBDAZ2-1.0D0-LAMBDA) + XP*XE*VIXEN5
     2          + XE**2*DLAMBDAY2*DLAMBDAZ2 - DPY2*DPZ2
         VIXEN8= YPP*(YE*VIXEN4-DMUZ2) + (YE*DMUY2-1.0D0-MU)*YE*DMUZ2 + ZPP*(ZE
     1          *VIXEN4-DMUY2) + ZE*DMUY2*(ZE*DMUZ2-1.0D0-MU) + XPP*XE*VIXEN4
     2          + XE**2*DMUY2*DMUZ2 - DQY2*DQZ2
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN7 + (VIXEN1/VIXEN2)*VIXEN8 + DPY2*DQZ2 + DQY2*DPZ2
         VIXEN6= YP*(YE*VIXEN4-DMUZ2) + (YE*DMUY2-1.0D0-MU)*YE*DLAMBDAZ2 + YPP*(YE*VIXEN5-DLAMBDAZ2)
     1          + (YE*DLAMBDAY2-1.0D0-LAMBDA)*YE*DMUZ2 + ZP*(ZE*VIXEN4-DMUY2) + ZE*DMUY2
     2          *(ZE*DLAMBDAZ2-1.0D0-LAMBDA) + ZPP*(ZE*VIXEN5-DLAMBDAY2) + ZE*DLAMBDAY2
     3          *(ZE*DMUZ2-1.0D0-MU) + XP*XE*VIXEN4 + XPP*XE*VIXEN5 + XE**2*(DMUY2*DLAMBDAZ2
     4          +DMUZ2*DLAMBDAY2)
         DNZ2= (V*VIXEN6 + DVZ2*DUY2 - U*VIXEN3 - DVY2*DUZ2)/V**2 - 2.0D0*DVZ2*N/V

         TORSHELL(3*B-1,3*B)= TORSHELL(3*B-1,3*B)+ M*DNZ2 + N*DMZ2
         TORSHELL(3*B,3*B-1)= TORSHELL(3*B,3*B-1)+ M*DNZ2 + N*DMZ2
C
C DX2DX3,DX3DX2
C
         N=(V*DUX2 - U*DVX2)/V**2
         VIXEN4= 1.0D0/RBC2 - (2.0D0*XE*(X(C)+X(D)-2.0D0*X(B))+2.0D0*((X(B)-X(C))*(X(B)-X(D))-MUCRAP))/RBC4 
     1          - 8.0D0*(X(B)-X(C))**2*MUCRAP/RBC6
         VIXEN5= 1.0D0/RBC2 - (2.0D0*XE*(X(C)+X(A)-2.0D0*X(B))+2.0D0*((X(B)-X(C))*(X(B)-X(A))-LAMBDACRAP))/RBC4 
     1          - 8.0D0*(X(B)-X(C))**2*LAMBDACRAP/RBC6
         VIXEN7= XP*(XE*VIXEN5+DLAMBDAX2-DLAMBDAX3) + (XE*DLAMBDAX2-1.0D0-LAMBDA)*(XE*DLAMBDAX3+LAMBDA)
     1          + (YP*YE+ZP*ZE)*VIXEN5 + (YE**2+ZE**2)*DLAMBDAX2*DLAMBDAX3 - DPX2*DPX3
         VIXEN8= XPP*(XE*VIXEN4+DMUX2-DMUX3) + (XE*DMUX2-1.0D0-MU)*(XE*DMUX3+MU)
     1          + (YPP*YE+ZPP*ZE)*VIXEN4 + (YE**2+ZE**2)*DMUX2*DMUX3 - DQX2*DQX3
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN7 + (VIXEN1/VIXEN2)*VIXEN8 + DPX2*DQX3 + DQX2*DPX3
         VIXEN6= XP*(XE*VIXEN4+DMUX2-DMUX3) + (XE*DMUX2-1.0D0-MU)*(XE*DLAMBDAX3+LAMBDA)
     1          + XPP*(XE*VIXEN5+DLAMBDAX2-DLAMBDAX3) + (XE*DLAMBDAX2-1.0D0-LAMBDA)*(XE*DMUX3+MU)
     2          + VIXEN4*(YP*YE+ZP*ZE) + VIXEN5*(YPP*YE+ZPP*ZE) 
     3          + (DMUX2*DLAMBDAX3+DMUX3*DLAMBDAX2)*(YE**2+ZE**2)
         DNX3= (V*VIXEN6 + DVX3*DUX2 - U*VIXEN3 - DUX3*DVX2)/V**2 - 2.0D0*DVX3*N/V

         TORSHELL(3*B-2,3*C-2)= TORSHELL(3*B-2,3*C-2)+ M*DNX3 + N*DMX3
         TORSHELL(3*C-2,3*B-2)= TORSHELL(3*C-2,3*B-2)+ M*DNX3 + N*DMX3
C
C DY2DY3,DY3DY2
C
         N=(V*DUY2 - U*DVY2)/V**2
         VIXEN4= 1.0D0/RBC2 - (2.0D0*YE*(Y(C)+Y(D)-2.0D0*Y(B))+2.0D0*((Y(B)-Y(C))*(Y(B)-Y(D))-MUCRAP))/RBC4 
     1          - 8.0D0*(Y(B)-Y(C))**2*MUCRAP/RBC6
         VIXEN5= 1.0D0/RBC2 - (2.0D0*YE*(Y(C)+Y(A)-2.0D0*Y(B))+2.0D0*((Y(B)-Y(C))*(Y(B)-Y(A))-LAMBDACRAP))/RBC4 
     1          - 8.0D0*(Y(B)-Y(C))**2*LAMBDACRAP/RBC6
         VIXEN7= YP*(YE*VIXEN5+DLAMBDAY2-DLAMBDAY3) + (YE*DLAMBDAY2-1.0D0-LAMBDA)*(YE*DLAMBDAY3+LAMBDA)
     1          + (XP*XE+ZP*ZE)*VIXEN5 + (XE**2+ZE**2)*DLAMBDAY2*DLAMBDAY3 - DPY2*DPY3
         VIXEN8= YPP*(YE*VIXEN4+DMUY2-DMUY3) + (YE*DMUY2-1.0D0-MU)*(YE*DMUY3+MU)
     1          + (XPP*XE+ZPP*ZE)*VIXEN4 + (XE**2+ZE**2)*DMUY2*DMUY3 - DQY2*DQY3
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN7 + (VIXEN1/VIXEN2)*VIXEN8 + DPY2*DQY3 + DQY2*DPY3
         VIXEN6= YP*(YE*VIXEN4+DMUY2-DMUY3) + (YE*DMUY2-1.0D0-MU)*(YE*DLAMBDAY3+LAMBDA)
     1          + YPP*(YE*VIXEN5+DLAMBDAY2-DLAMBDAY3) + (YE*DLAMBDAY2-1.0D0-LAMBDA)*(YE*DMUY3+MU)
     2          + VIXEN4*(XP*XE+ZP*ZE) + VIXEN5*(XPP*XE+ZPP*ZE) 
     3          + (DMUY2*DLAMBDAY3+DMUY3*DLAMBDAY2)*(XE**2+ZE**2)
         DNY3= (V*VIXEN6 + DVY3*DUY2 - U*VIXEN3 - DUY3*DVY2)/V**2 - 2.0D0*DVY3*N/V

         TORSHELL(3*B-1,3*C-1)= TORSHELL(3*B-1,3*C-1)+ M*DNY3 + N*DMY3
         TORSHELL(3*C-1,3*B-1)= TORSHELL(3*C-1,3*B-1)+ M*DNY3 + N*DMY3
C
C DZ2DZ3,DZ3DZ2
C
         N=(V*DUZ2 - U*DVZ2)/V**2
         VIXEN4= 1.0D0/RBC2 - (2.0D0*ZE*(Z(C)+Z(D)-2.0D0*Z(B))+2.0D0*((Z(B)-Z(C))*(Z(B)-Z(D))-MUCRAP))/RBC4 
     1          - 8.0D0*(Z(B)-Z(C))**2*MUCRAP/RBC6
         VIXEN5= 1.0D0/RBC2 - (2.0D0*ZE*(Z(C)+Z(A)-2.0D0*Z(B))+2.0D0*((Z(B)-Z(C))*(Z(B)-Z(A))-LAMBDACRAP))/RBC4 
     1          - 8.0D0*(Z(B)-Z(C))**2*LAMBDACRAP/RBC6
         VIXEN7= ZP*(ZE*VIXEN5+DLAMBDAZ2-DLAMBDAZ3) + (ZE*DLAMBDAZ2-1.0D0-LAMBDA)*(ZE*DLAMBDAZ3+LAMBDA)
     1          + (XP*XE+YP*YE)*VIXEN5 + (XE**2+YE**2)*DLAMBDAZ2*DLAMBDAZ3 - DPZ2*DPZ3
         VIXEN8= ZPP*(ZE*VIXEN4+DMUZ2-DMUZ3) + (ZE*DMUZ2-1.0D0-MU)*(ZE*DMUZ3+MU)
     1          + (XPP*XE+YPP*YE)*VIXEN4 + (XE**2+YE**2)*DMUZ2*DMUZ3 - DQZ2*DQZ3
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN7 + (VIXEN1/VIXEN2)*VIXEN8 + DPZ2*DQZ3 + DQZ2*DPZ3
         VIXEN6= ZP*(ZE*VIXEN4+DMUZ2-DMUZ3) + (ZE*DMUZ2-1.0D0-MU)*(ZE*DLAMBDAZ3+LAMBDA)
     1          + ZPP*(ZE*VIXEN5+DLAMBDAZ2-DLAMBDAZ3) + (ZE*DLAMBDAZ2-1-LAMBDA)*(ZE*DMUZ3+MU)
     2          + VIXEN4*(XP*XE+YP*YE) + VIXEN5*(XPP*XE+YPP*YE) 
     3          + (DMUZ2*DLAMBDAZ3+DMUZ3*DLAMBDAZ2)*(XE**2+YE**2)
         DNZ3= (V*VIXEN6 + DVZ3*DUZ2 - U*VIXEN3 - DUZ3*DVZ2)/V**2 - 2.0D0*DVZ3*N/V

         TORSHELL(3*B,3*C)= TORSHELL(3*B,3*C)+ M*DNZ3 + N*DMZ3
         TORSHELL(3*C,3*B)= TORSHELL(3*C,3*B)+ M*DNZ3 + N*DMZ3
C
C DX2DY3,DY3DX2
C
         N=(V*DUX2 - U*DVX2)/V**2
         VIXEN4= 8.0D0*(X(B)-X(C))*YE*MUCRAP/RBC6 - 2.0D0*((X(C)+X(D)-2*X(B))*YE
     1          +(X(B)-X(C))*(Y(B)-Y(D)))/RBC4
         VIXEN5= 8.0D0*(X(B)-X(C))*YE*LAMBDACRAP/RBC6 - 2.0D0*((X(C)+X(A)-2*X(B))*YE
     1          +(X(B)-X(C))*(Y(B)-Y(A)))/RBC4
         VIXEN7= XP*(XE*VIXEN5-DLAMBDAY3) + (XE*DLAMBDAX2-1.0D0-LAMBDA)*XE*DLAMBDAY3 + 
     1          YP*(YE*VIXEN5+DLAMBDAX2) + YE*DLAMBDAX2*(YE*DLAMBDAY3+LAMBDA) 
     2          + ZP*ZE*VIXEN5 + ZE**2*DLAMBDAX2*DLAMBDAY3 - DPX2*DPY3
         VIXEN8= XPP*(XE*VIXEN4-DMUY3) + (XE*DMUX2-1.0D0-MU)*XE*DMUY3 + 
     1          YPP*(YE*VIXEN4+DMUX2) + YE*DMUX2*(YE*DMUY3+MU) + 
     2          ZPP*ZE*VIXEN4 + ZE**2*DMUX2*DMUY3 - DQX2*DQY3
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN7 + (VIXEN1/VIXEN2)*VIXEN8 + DPX2*DQY3 + DQX2*DPY3
         VIXEN6= XP*(XE*VIXEN4-DMUY3) + (XE*DMUX2-1.0D0-MU)*XE*DLAMBDAY3
     1          + XPP*(XE*VIXEN5-DLAMBDAY3) + (XE*DLAMBDAX2-1.0D0-LAMBDA)*XE*DMUY3
     2          + YP*(YE*VIXEN4+DMUX2) + YE*DMUX2*(YE*DLAMBDAY3+LAMBDA)
     3          + YPP*(YE*VIXEN5+DLAMBDAX2) + YE*DLAMBDAX2*(YE*DMUY3+MU)
     4          + ZP*ZE*VIXEN4 + ZE**2*(DMUX2*DLAMBDAY3+DLAMBDAX2*DMUY3) + ZPP*ZE*VIXEN5
         DNY3= (V*VIXEN6 + DVY3*DUX2 - U*VIXEN3 - DUY3*DVX2)/V**2 - 2.0D0*DVY3*N/V

         TORSHELL(3*B-2,3*C-1)= TORSHELL(3*B-2,3*C-1)+ M*DNY3 + N*DMY3
         TORSHELL(3*C-1,3*B-2)= TORSHELL(3*C-1,3*B-2)+ M*DNY3 + N*DMY3
C
C DX2DZ3,DZ3DX2
C
         N=(V*DUX2 - U*DVX2)/V**2
         VIXEN4= 8.0D0*(X(B)-X(C))*ZE*MUCRAP/RBC6 - 2.0D0*((X(C)+X(D)-2*X(B))*ZE
     1          +(X(B)-X(C))*(Z(B)-Z(D)))/RBC4
         VIXEN5= 8.0D0*(X(B)-X(C))*ZE*LAMBDACRAP/RBC6 - 2.0D0*((X(C)+X(A)-2*X(B))*ZE
     1          +(X(B)-X(C))*(Z(B)-Z(A)))/RBC4
         VIXEN7= XP*(XE*VIXEN5-DLAMBDAZ3) + (XE*DLAMBDAX2-1.0D0-LAMBDA)*XE*DLAMBDAZ3 + 
     1          ZP*(ZE*VIXEN5+DLAMBDAX2) + ZE*DLAMBDAX2*(ZE*DLAMBDAZ3+LAMBDA) 
     2          + YP*YE*VIXEN5 + YE**2*DLAMBDAX2*DLAMBDAZ3 - DPX2*DPZ3
         VIXEN8= XPP*(XE*VIXEN4-DMUZ3) + (XE*DMUX2-1.0D0-MU)*XE*DMUZ3 + 
     1          ZPP*(ZE*VIXEN4+DMUX2) + ZE*DMUX2*(ZE*DMUZ3+MU) + 
     2          YPP*YE*VIXEN4 + YE**2*DMUX2*DMUZ3 - DQX2*DQZ3
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN7 + (VIXEN1/VIXEN2)*VIXEN8 + DPX2*DQZ3 + DQX2*DPZ3
         VIXEN6= XP*(XE*VIXEN4-DMUZ3) + (XE*DMUX2-1.0D0-MU)*XE*DLAMBDAZ3
     1          + XPP*(XE*VIXEN5-DLAMBDAZ3) + (XE*DLAMBDAX2-1.0D0-LAMBDA)*XE*DMUZ3
     2          + ZP*(ZE*VIXEN4+DMUX2) + ZE*DMUX2*(ZE*DLAMBDAZ3+LAMBDA)
     3          + ZPP*(ZE*VIXEN5+DLAMBDAX2) + ZE*DLAMBDAX2*(ZE*DMUZ3+MU)
     4          + YP*YE*VIXEN4 + YE**2*(DMUX2*DLAMBDAZ3+DLAMBDAX2*DMUZ3) + YPP*YE*VIXEN5
         DNZ3= (V*VIXEN6 + DVZ3*DUX2 - U*VIXEN3 - DUZ3*DVX2)/V**2 - 2.0D0*DVZ3*N/V

         TORSHELL(3*B-2,3*C)= TORSHELL(3*B-2,3*C)+ M*DNZ3 + N*DMZ3
         TORSHELL(3*C,3*B-2)= TORSHELL(3*C,3*B-2)+ M*DNZ3 + N*DMZ3
C
C DY2DZ3,DZ3DY2
C
         N=(V*DUY2 - U*DVY2)/V**2
         VIXEN4= 8.0D0*(Y(B)-Y(C))*ZE*MUCRAP/RBC6 - 2.0D0*((Y(C)+Y(D)-2*Y(B))*ZE
     1          +(Y(B)-Y(C))*(Z(B)-Z(D)))/RBC4
         VIXEN5= 8.0D0*(Y(B)-Y(C))*ZE*LAMBDACRAP/RBC6 - 2.0D0*((Y(C)+Y(A)-2*Y(B))*ZE
     1          +(Y(B)-Y(C))*(Z(B)-Z(A)))/RBC4
         VIXEN7= YP*(YE*VIXEN5-DLAMBDAZ3) + (YE*DLAMBDAY2-1.0D0-LAMBDA)*YE*DLAMBDAZ3 +
     1          ZP*(ZE*VIXEN5+DLAMBDAY2) + ZE*DLAMBDAY2*(ZE*DLAMBDAZ3+LAMBDA)
     2          + XP*XE*VIXEN5 + XE**2*DLAMBDAY2*DLAMBDAZ3 - DPY2*DPZ3
         VIXEN8= YPP*(YE*VIXEN4-DMUZ3) + (YE*DMUY2-1.0D0-MU)*YE*DMUZ3 +
     1          ZPP*(ZE*VIXEN4+DMUY2) + ZE*DMUY2*(ZE*DMUZ3+MU) +
     2          XPP*XE*VIXEN4 + XE**2*DMUY2*DMUZ3 - DQY2*DQZ3
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN7 + (VIXEN1/VIXEN2)*VIXEN8 + DPY2*DQZ3 + DQY2*DPZ3
         VIXEN6= YP*(YE*VIXEN4-DMUZ3) + (YE*DMUY2-1.0D0-MU)*YE*DLAMBDAZ3
     1          + YPP*(YE*VIXEN5-DLAMBDAZ3) + (YE*DLAMBDAY2-1.0D0-LAMBDA)*YE*DMUZ3
     2          + ZP*(ZE*VIXEN4+DMUY2) + ZE*DMUY2*(ZE*DLAMBDAZ3+LAMBDA)
     3          + ZPP*(ZE*VIXEN5+DLAMBDAY2) + ZE*DLAMBDAY2*(ZE*DMUZ3+MU)
     4          + XP*XE*VIXEN4 + XE**2*(DMUY2*DLAMBDAZ3+DLAMBDAY2*DMUZ3) + XPP*XE*VIXEN5
         DNZ3= (V*VIXEN6 + DVZ3*DUY2 - U*VIXEN3 - DUZ3*DVY2)/V**2 - 2.0D0*DVZ3*N/V

         TORSHELL(3*B-1,3*C)= TORSHELL(3*B-1,3*C)+ M*DNZ3 + N*DMZ3
         TORSHELL(3*C,3*B-1)= TORSHELL(3*C,3*B-1)+ M*DNZ3 + N*DMZ3
C
C DY2DX3,DX3DY2
C
         N=(V*DUY2 - U*DVY2)/V**2
         VIXEN4= 8.0D0*(Y(B)-Y(C))*XE*MUCRAP/RBC6 - 2.0D0*((Y(C)+Y(D)-2*Y(B))*XE
     1          +(Y(B)-Y(C))*(X(B)-X(D)))/RBC4
         VIXEN5= 8.0D0*(Y(B)-Y(C))*XE*LAMBDACRAP/RBC6 - 2.0D0*((Y(C)+Y(A)-2*Y(B))*XE
     1          +(Y(B)-Y(C))*(X(B)-X(A)))/RBC4
         VIXEN7= YP*(YE*VIXEN5-DLAMBDAX3) + (YE*DLAMBDAY2-1.0D0-LAMBDA)*YE*DLAMBDAX3 +
     1          XP*(XE*VIXEN5+DLAMBDAY2) + XE*DLAMBDAY2*(XE*DLAMBDAX3+LAMBDA)
     2          + ZP*ZE*VIXEN5 + ZE**2*DLAMBDAY2*DLAMBDAX3 - DPY2*DPX3
         VIXEN8= YPP*(YE*VIXEN4-DMUX3) + (YE*DMUY2-1.0D0-MU)*YE*DMUX3 +
     1          XPP*(XE*VIXEN4+DMUY2) + XE*DMUY2*(XE*DMUX3+MU) +
     2          ZPP*ZE*VIXEN4 + ZE**2*DMUY2*DMUX3 - DQY2*DQX3
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN7 + (VIXEN1/VIXEN2)*VIXEN8 + DPY2*DQX3 + DQY2*DPX3
         VIXEN6= YP*(YE*VIXEN4-DMUX3) + (YE*DMUY2-1.0D0-MU)*YE*DLAMBDAX3
     1          + YPP*(YE*VIXEN5-DLAMBDAX3) + (YE*DLAMBDAY2-1.0D0-LAMBDA)*YE*DMUX3
     2          + XP*(XE*VIXEN4+DMUY2) + XE*DMUY2*(XE*DLAMBDAX3+LAMBDA)
     3          + XPP*(XE*VIXEN5+DLAMBDAY2) + XE*DLAMBDAY2*(XE*DMUX3+MU)
     4          + ZP*ZE*VIXEN4 + ZE**2*(DMUY2*DLAMBDAX3+DLAMBDAY2*DMUX3) + ZPP*ZE*VIXEN5
         DNX3= (V*VIXEN6 + DVX3*DUY2 - U*VIXEN3 - DUX3*DVY2)/V**2 - 2.0D0*DVX3*N/V

         TORSHELL(3*B-1,3*C-2)= TORSHELL(3*B-1,3*C-2)+ M*DNX3 + N*DMX3
         TORSHELL(3*C-2,3*B-1)= TORSHELL(3*C-2,3*B-1)+ M*DNX3 + N*DMX3
C
C DZ2DX3,DX3DZ2
C
         N=(V*DUZ2 - U*DVZ2)/V**2
         VIXEN4= 8.0D0*(Z(B)-Z(C))*XE*MUCRAP/RBC6 - 2.0D0*((Z(C)+Z(D)-2*Z(B))*XE
     1          +(Z(B)-Z(C))*(X(B)-X(D)))/RBC4
         VIXEN5= 8.0D0*(Z(B)-Z(C))*XE*LAMBDACRAP/RBC6 - 2.0D0*((Z(C)+Z(A)-2*Z(B))*XE
     1          +(Z(B)-Z(C))*(X(B)-X(A)))/RBC4
         VIXEN7= ZP*(ZE*VIXEN5-DLAMBDAX3) + (ZE*DLAMBDAZ2-1.0D0-LAMBDA)*ZE*DLAMBDAX3 +
     1          XP*(XE*VIXEN5+DLAMBDAZ2) + XE*DLAMBDAZ2*(XE*DLAMBDAX3+LAMBDA)
     2          + YP*YE*VIXEN5 + YE**2*DLAMBDAZ2*DLAMBDAX3 - DPZ2*DPX3
         VIXEN8= ZPP*(ZE*VIXEN4-DMUX3) + (ZE*DMUZ2-1.0D0-MU)*ZE*DMUX3 +
     1          XPP*(XE*VIXEN4+DMUZ2) + XE*DMUZ2*(XE*DMUX3+MU) +
     2          YPP*YE*VIXEN4 + YE**2*DMUZ2*DMUX3 - DQZ2*DQX3
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN7 + (VIXEN1/VIXEN2)*VIXEN8 + DPZ2*DQX3 + DQZ2*DPX3
         VIXEN6= ZP*(ZE*VIXEN4-DMUX3) + (ZE*DMUZ2-1.0D0-MU)*ZE*DLAMBDAX3
     1          + ZPP*(ZE*VIXEN5-DLAMBDAX3) + (ZE*DLAMBDAZ2-1.0D0-LAMBDA)*ZE*DMUX3
     2          + XP*(XE*VIXEN4+DMUZ2) + XE*DMUZ2*(XE*DLAMBDAX3+LAMBDA)
     3          + XPP*(XE*VIXEN5+DLAMBDAZ2) + XE*DLAMBDAZ2*(XE*DMUX3+MU)
     4          + YP*YE*VIXEN4 + YE**2*(DMUZ2*DLAMBDAX3+DLAMBDAZ2*DMUX3) + YPP*YE*VIXEN5
         DNX3= (V*VIXEN6 + DVX3*DUZ2 - U*VIXEN3 - DUX3*DVZ2)/V**2 - 2.0D0*DVX3*N/V

         TORSHELL(3*B,3*C-2)= TORSHELL(3*B,3*C-2)+ M*DNX3 + N*DMX3
         TORSHELL(3*C-2,3*B)= TORSHELL(3*C-2,3*B)+ M*DNX3 + N*DMX3
C
C DZ2DY3,DY3DZ2
C
         N=(V*DUZ2 - U*DVZ2)/V**2
         VIXEN4= 8.0D0*(Z(B)-Z(C))*YE*MUCRAP/RBC6 - 2.0D0*((Z(C)+Z(D)-2*Z(B))*YE
     1          +(Z(B)-Z(C))*(Y(B)-Y(D)))/RBC4
         VIXEN5= 8.0D0*(Z(B)-Z(C))*YE*LAMBDACRAP/RBC6 - 2.0D0*((Z(C)+Z(A)-2*Z(B))*YE
     1          +(Z(B)-Z(C))*(Y(B)-Y(A)))/RBC4
         VIXEN7= ZP*(ZE*VIXEN5-DLAMBDAY3) + (ZE*DLAMBDAZ2-1.0D0-LAMBDA)*ZE*DLAMBDAY3 +
     1          YP*(YE*VIXEN5+DLAMBDAZ2) + YE*DLAMBDAZ2*(YE*DLAMBDAY3+LAMBDA)
     2          + XP*XE*VIXEN5 + XE**2*DLAMBDAZ2*DLAMBDAY3 - DPZ2*DPY3
         VIXEN8= ZPP*(ZE*VIXEN4-DMUY3) + (ZE*DMUZ2-1.0D0-MU)*ZE*DMUY3 +
     1          YPP*(YE*VIXEN4+DMUZ2) + YE*DMUZ2*(YE*DMUY3+MU) +
     2          XPP*XE*VIXEN4 + XE**2*DMUZ2*DMUY3 - DQZ2*DQY3
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN7 + (VIXEN1/VIXEN2)*VIXEN8 + DPZ2*DQY3 + DQZ2*DPY3
         VIXEN6= ZP*(ZE*VIXEN4-DMUY3) + (ZE*DMUZ2-1.0D0-MU)*ZE*DLAMBDAY3
     1          + ZPP*(ZE*VIXEN5-DLAMBDAY3) + (ZE*DLAMBDAZ2-1.0D0-LAMBDA)*ZE*DMUY3
     2          + YP*(YE*VIXEN4+DMUZ2) + YE*DMUZ2*(YE*DLAMBDAY3+LAMBDA)
     3          + YPP*(YE*VIXEN5+DLAMBDAZ2) + YE*DLAMBDAZ2*(YE*DMUY3+MU)
     4          + XP*XE*VIXEN4 + XE**2*(DMUZ2*DLAMBDAY3+DLAMBDAZ2*DMUY3) + XPP*XE*VIXEN5
         DNY3= (V*VIXEN6 + DVY3*DUZ2 - U*VIXEN3 - DUY3*DVZ2)/V**2 - 2.0D0*DVY3*N/V

         TORSHELL(3*B,3*C-1)= TORSHELL(3*B,3*C-1)+ M*DNY3 + N*DMY3
         TORSHELL(3*C-1,3*B)= TORSHELL(3*C-1,3*B)+ M*DNY3 + N*DMY3



      END DO

C      PRINT *,' '
C      PRINT *,'TORSION TEST....'
C      PRINT *,' '
C      DO I=1,ATOMS
C         VIXEN1=DTORSEBYDX(I)-TORSCRAP(3*I-2)
C         WRITE(*,FMT='(A,I2,A,F20.10,5X,F20.10,A,F20.10)') 'X(',I,')',DTORSEBYDX(I),TORSCRAP(3*I-2),' DIFF ',VIXEN1
C         VIXEN1=DTORSEBYDY(I)-TORSCRAP(3*I-1)
C         WRITE(*,FMT='(A,I2,A,F20.10,5X,F20.10,A,F20.10)') 'Y(',I,')',DTORSEBYDY(I),TORSCRAP(3*I-1),' DIFF ',VIXEN1
C         VIXEN1=DTORSEBYDZ(I)-TORSCRAP(3*I)
C         WRITE(*,FMT='(A,I2,A,F20.10,5X,F20.10,A,F20.10)') 'Z(',I,')',DTORSEBYDZ(I),TORSCRAP(3*I),' DIFF ',VIXEN1
C      END DO
C
C************************************
C IMPROPER TORSION ANGLE DERIVATIVES*
C************************************
C
      DO I=1,IMP
         A=IA1(I)
         B=IA2(I)
         C=IA3(I)
         D=IA4(I)
C         PRINT *,A,B,C,D,' : ',TYPECH(A),'-',TYPECH(B),'-',TYPECH(C),'-',TYPECH(D)

         RBC=R(B,C)
         RBC2=RBC**2
         RBC4=RBC**4
         RBC6=RBC**6
         XE=X(C)-X(B)
         YE=Y(C)-Y(B)
         ZE=Z(C)-Z(B)
         LAMBDA=((X(B)-X(A))*XE+(Y(B)-Y(A))*YE+(Z(B)-Z(A))*ZE)/RBC2
         MU=    ((X(B)-X(D))*XE+(Y(B)-Y(D))*YE+(Z(B)-Z(D))*ZE)/RBC2
         XP= X(A)+LAMBDA*X(C)-(1.0D0+LAMBDA)*X(B)
         YP= Y(A)+LAMBDA*Y(C)-(1.0D0+LAMBDA)*Y(B)
         ZP= Z(A)+LAMBDA*Z(C)-(1.0D0+LAMBDA)*Z(B)
         XPP= X(D)+MU*X(C)-(1.0D0+MU)*X(B)
         YPP= Y(D)+MU*Y(C)-(1.0D0+MU)*Y(B)
         ZPP= Z(D)+MU*Z(C)-(1.0D0+MU)*Z(B)
         U= XP*XPP + YP*YPP + ZP*ZPP
         V= SQRT(XP**2 + YP**2 + ZP**2)*SQRT(XPP**2 + YPP**2 + ZPP**2)

         LAMBDACRAP= (X(B)-X(A))*XE + (Y(B)-Y(A))*YE + (Z(B)-Z(A))*ZE
         MUCRAP= (X(B)-X(D))*XE + (Y(B)-Y(D))*YE + (Z(B)-Z(D))*ZE

         DLAMBDAX1= (X(B)-X(C))/RBC2
         DLAMBDAY1= (Y(B)-Y(C))/RBC2
         DLAMBDAZ1= (Z(B)-Z(C))/RBC2
         DLAMBDAX2= (X(C)+X(A)-2.0D0*X(B))/RBC2 - 2.0D0*(X(B)-X(C))*LAMBDACRAP/RBC4 
         DLAMBDAY2= (Y(C)+Y(A)-2.0D0*Y(B))/RBC2 - 2.0D0*(Y(B)-Y(C))*LAMBDACRAP/RBC4
         DLAMBDAZ2= (Z(C)+Z(A)-2.0D0*Z(B))/RBC2 - 2.0D0*(Z(B)-Z(C))*LAMBDACRAP/RBC4
         DLAMBDAX3= (X(B)-X(A))/RBC2 - 2.0D0*XE*LAMBDACRAP/RBC4
         DLAMBDAY3= (Y(B)-Y(A))/RBC2 - 2.0D0*YE*LAMBDACRAP/RBC4
         DLAMBDAZ3= (Z(B)-Z(A))/RBC2 - 2.0D0*ZE*LAMBDACRAP/RBC4
         DMUX2= (X(C)+X(D)-2.0D0*X(B))/RBC2 - 2.0D0*(X(B)-X(C))*MUCRAP/RBC4
         DMUY2= (Y(C)+Y(D)-2.0D0*Y(B))/RBC2 - 2.0D0*(Y(B)-Y(C))*MUCRAP/RBC4
         DMUZ2= (Z(C)+Z(D)-2.0D0*Z(B))/RBC2 - 2.0D0*(Z(B)-Z(C))*MUCRAP/RBC4
         DMUX3= (X(B)-X(D))/RBC2 - 2.0D0*XE*MUCRAP/RBC4
         DMUY3= (Y(B)-Y(D))/RBC2 - 2.0D0*YE*MUCRAP/RBC4
         DMUZ3= (Z(B)-Z(D))/RBC2 - 2.0D0*ZE*MUCRAP/RBC4
         DMUX4= (X(B)-X(C))/RBC2
         DMUY4= (Y(B)-Y(C))/RBC2
         DMUZ4= (Z(B)-Z(C))/RBC2

         VIXEN1=SQRT(XP**2 + YP**2 + ZP**2)
         VIXEN2=SQRT(XPP**2 + YPP**2 + ZPP**2)
         DPX1=(1.0D0/VIXEN1)*(XP*(1.0D0+XE*DLAMBDAX1) + YP*YE*DLAMBDAX1 + ZP*ZE*DLAMBDAX1)
         DVX1=VIXEN2*DPX1
         DUX1=XPP*(1.0D0+XE*DLAMBDAX1) + YPP*YE*DLAMBDAX1 + ZPP*ZE*DLAMBDAX1

         DPY1=(1.0D0/VIXEN1)*(YP*(1.0D0+YE*DLAMBDAY1) + XP*XE*DLAMBDAY1 + ZP*ZE*DLAMBDAY1)
         DVY1=VIXEN2*DPY1
         DUY1=YPP*(1.0D0+YE*DLAMBDAY1) + XPP*XE*DLAMBDAY1 + ZPP*ZE*DLAMBDAY1

         DPZ1=(1.0D0/VIXEN1)*(ZP*(1.0D0+ZE*DLAMBDAZ1) + XP*XE*DLAMBDAZ1 + YP*YE*DLAMBDAZ1)
         DVZ1=VIXEN2*DPZ1
         DUZ1=ZPP*(1.0D0+ZE*DLAMBDAZ1) + XPP*XE*DLAMBDAZ1 + YPP*YE*DLAMBDAZ1

         DUX2=XP*(XE*DMUX2-1.0D0-MU) + XPP*(XE*DLAMBDAX2-1.0D0-LAMBDA) + YE*(YP*DMUX2 + YPP*DLAMBDAX2) +
     1          ZE*(ZP*DMUX2 + ZPP*DLAMBDAX2)
         DPX2=(1.0D0/VIXEN1)*(XP*(XE*DLAMBDAX2-1.0D0-LAMBDA) + YP*YE*DLAMBDAX2 + ZP*ZE*DLAMBDAX2)
         VIXEN6=VIXEN2*DPX2
         DQX2=(1.0D0/VIXEN2)*(XPP*(XE*DMUX2-1.0D0-MU) + YPP*YE*DMUX2 + ZPP*ZE*DMUX2)
         VIXEN7=VIXEN1*DQX2
         DVX2=VIXEN6+VIXEN7

         DUY2=YP*(YE*DMUY2-1.0D0-MU) + YPP*(YE*DLAMBDAY2-1.0D0-LAMBDA) + XE*(XP*DMUY2 + XPP*DLAMBDAY2) +
     1          ZE*(ZP*DMUY2 + ZPP*DLAMBDAY2)
         DPY2=(1.0D0/VIXEN1)*(YP*(YE*DLAMBDAY2-1.0D0-LAMBDA) + XP*XE*DLAMBDAY2 + ZP*ZE*DLAMBDAY2)
         VIXEN6=VIXEN2*DPY2
         DQY2=(1.0D0/VIXEN2)*(YPP*(YE*DMUY2-1.0D0-MU) + XPP*XE*DMUY2 + ZPP*ZE*DMUY2)
         VIXEN7=VIXEN1*DQY2
         DVY2=VIXEN6+VIXEN7

         DUZ2=ZP*(ZE*DMUZ2-1.0D0-MU) + ZPP*(ZE*DLAMBDAZ2-1.0D0-LAMBDA) + XE*(XP*DMUZ2 + XPP*DLAMBDAZ2) +
     1          YE*(YP*DMUZ2 + YPP*DLAMBDAZ2)
         DPZ2=(1.0D0/VIXEN1)*(ZP*(ZE*DLAMBDAZ2-1.0D0-LAMBDA) + XP*XE*DLAMBDAZ2 + YP*YE*DLAMBDAZ2)
         VIXEN6=VIXEN2*DPZ2
         DQZ2=(1.0D0/VIXEN2)*(ZPP*(ZE*DMUZ2-1.0D0-MU) + XPP*XE*DMUZ2 + YPP*YE*DMUZ2)
         VIXEN7=VIXEN1*DQZ2
         DVZ2=VIXEN6+VIXEN7

         DUX3=XP*(MU + XE*DMUX3) + XPP*(LAMBDA + XE*DLAMBDAX3) + YE*(YP*DMUX3 + YPP*DLAMBDAX3) +
     1          ZE*(ZP*DMUX3 + ZPP*DLAMBDAX3)
         DPX3=(1.0D0/VIXEN1)*(XP*(LAMBDA + XE*DLAMBDAX3) + YP*YE*DLAMBDAX3 + ZP*ZE*DLAMBDAX3)
         VIXEN6=VIXEN2*DPX3
         DQX3=(1.0D0/VIXEN2)*(XPP*(MU + XE*DMUX3) + YPP*YE*DMUX3 + ZPP*ZE*DMUX3)
         VIXEN7=VIXEN1*DQX3
         DVX3=VIXEN6+VIXEN7

         DUY3=YP*(MU + YE*DMUY3) + YPP*(LAMBDA + YE*DLAMBDAY3) + XE*(XP*DMUY3 + XPP*DLAMBDAY3) +
     1          ZE*(ZP*DMUY3 + ZPP*DLAMBDAY3)
         DPY3=(1.0D0/VIXEN1)*(YP*(LAMBDA + YE*DLAMBDAY3) + XP*XE*DLAMBDAY3 + ZP*ZE*DLAMBDAY3)
         VIXEN6=VIXEN2*DPY3
         DQY3=(1.0D0/VIXEN2)*(YPP*(MU + YE*DMUY3) + XPP*XE*DMUY3 + ZPP*ZE*DMUY3)
         VIXEN7=VIXEN1*DQY3
         DVY3=VIXEN6+VIXEN7

         DUZ3=ZP*(MU + ZE*DMUZ3) + ZPP*(LAMBDA + ZE*DLAMBDAZ3) + XE*(XP*DMUZ3 + XPP*DLAMBDAZ3) +
     1          YE*(YP*DMUZ3 + YPP*DLAMBDAZ3)
         DPZ3=(1.0D0/VIXEN1)*(ZP*(LAMBDA + ZE*DLAMBDAZ3) + XP*XE*DLAMBDAZ3 + YP*YE*DLAMBDAZ3)
         VIXEN6=VIXEN2*DPZ3
         DQZ3=(1.0D0/VIXEN2)*(ZPP*(MU + ZE*DMUZ3) + XPP*XE*DMUZ3 + YPP*YE*DMUZ3)
         VIXEN7=VIXEN1*DQZ3
         DVZ3=VIXEN6+VIXEN7

         DUX4=XP*(1.0D0 + XE*DMUX4) + YP*YE*DMUX4 + ZP*ZE*DMUX4
         DQX4=(1.0D0/VIXEN2)*(XPP*(1.0D0 + XE*DMUX4) + YPP*YE*DMUX4 + ZPP*ZE*DMUX4)
         DVX4=VIXEN1*DQX4
  
         DUY4=YP*(1.0D0 + YE*DMUY4) + XP*XE*DMUY4 + ZP*ZE*DMUY4
         DQY4=(1.0D0/VIXEN2)*(YPP*(1.0D0 + YE*DMUY4) + XPP*XE*DMUY4 + ZPP*ZE*DMUY4)
         DVY4=VIXEN1*DQY4

         DUZ4=ZP*(1.0D0 + ZE*DMUZ4) + XP*XE*DMUZ4 + YP*YE*DMUZ4
         DQZ4=(1.0D0/VIXEN2)*(ZPP*(1.0D0 + ZE*DMUZ4) + XPP*XE*DMUZ4 + YPP*YE*DMUZ4)
         DVZ4=VIXEN1*DQZ4
C
C NOW SECOND DERIVATIVES
C
         CALL HAIRYIMPHELL(U,V)
C
C DX1DX1
C
         N=(V*DUX1 - U*DVX1)/V**2
         VIXEN3= ((1.0D0+XE*DLAMBDAX1)**2 + YE**2*DLAMBDAX1**2 + ZE**2*DLAMBDAX1**2 - DPX1**2)/VIXEN1
         DNX1= -(U*VIXEN2/V**2)*VIXEN3 -2.0D0*DVX1*(V*DUX1 - U*DVX1)/V**3

         IMPHELL(3*A-2,3*A-2) =IMPHELL(3*A-2,3*A-2) + M*DNX1 + N*DMX1
C
C DY1DY1
C
         N=(V*DUY1 - U*DVY1)/V**2
         VIXEN3= ((1.0D0+YE*DLAMBDAY1)**2 + XE**2*DLAMBDAY1**2 + ZE**2*DLAMBDAY1**2 - DPY1**2)/VIXEN1
         DNY1= -(U*VIXEN2/V**2)*VIXEN3 -2.0D0*DVY1*(V*DUY1 - U*DVY1)/V**3

         IMPHELL(3*A-1,3*A-1) =IMPHELL(3*A-1,3*A-1) + M*DNY1 + N*DMY1
C
C DZ1DZ1
C
         N=(V*DUZ1 - U*DVZ1)/V**2
         VIXEN3= ((1.0D0+ZE*DLAMBDAZ1)**2 + XE**2*DLAMBDAZ1**2 + YE**2*DLAMBDAZ1**2 - DPZ1**2)/VIXEN1
         DNZ1= -(U*VIXEN2/V**2)*VIXEN3 -2.0D0*DVZ1*(V*DUZ1 - U*DVZ1)/V**3

         IMPHELL(3*A,3*A) =IMPHELL(3*A,3*A) + M*DNZ1 + N*DMZ1
C
C DX1DY1,DY1DX1
C
         N=(V*DUX1 - U*DVX1)/V**2
         VIXEN3= ((1.0D0+XE*DLAMBDAX1)*XE*DLAMBDAY1 + YE*DLAMBDAX1*(1.0D0+YE*DLAMBDAY1) 
     1          +ZE**2*DLAMBDAX1*DLAMBDAY1 - DPX1*DPY1)/VIXEN1
         DNY1= (DUX1*DVY1 - DUY1*DVX1 - U*VIXEN2*VIXEN3)/V**2 - 2.0D0*DVY1*(V*DUX1 - U*DVX1)/V**3
         
         IMPHELL(3*A-2,3*A-1)= IMPHELL(3*A-2,3*A-1) + M*DNY1 + N*DMY1
         IMPHELL(3*A-1,3*A-2)= IMPHELL(3*A-1,3*A-2) + M*DNY1 + N*DMY1
C
C DX1DZ1,DZ1DX1
C
         N=(V*DUX1 - U*DVX1)/V**2
         VIXEN3= ((1.0D0+XE*DLAMBDAX1)*XE*DLAMBDAZ1 + ZE*DLAMBDAX1*(1.0D0+ZE*DLAMBDAZ1)
     1          +YE**2*DLAMBDAX1*DLAMBDAZ1 - DPX1*DPZ1)/VIXEN1
         DNZ1= (DUX1*DVZ1 - DUZ1*DVX1 - U*VIXEN2*VIXEN3)/V**2 - 2.0D0*DVZ1*(V*DUX1 - U*DVX1)/V**3

         IMPHELL(3*A-2,3*A)= IMPHELL(3*A-2,3*A) + M*DNZ1 + N*DMZ1
         IMPHELL(3*A,3*A-2)= IMPHELL(3*A,3*A-2) + M*DNZ1 + N*DMZ1
C
C DY1DZ1,DZ1DY1
C
         N=(V*DUY1 - U*DVY1)/V**2
         VIXEN3= ((1.0D0+YE*DLAMBDAY1)*YE*DLAMBDAZ1 + ZE*DLAMBDAY1*(1.0D0+ZE*DLAMBDAZ1)
     1          +XE**2*DLAMBDAY1*DLAMBDAZ1 - DPY1*DPZ1)/VIXEN1
         DNZ1= (DUY1*DVZ1 - DUZ1*DVY1 - U*VIXEN2*VIXEN3)/V**2 - 2.0D0*DVZ1*(V*DUY1 - U*DVY1)/V**3

         IMPHELL(3*A-1,3*A)= IMPHELL(3*A-1,3*A) + M*DNZ1 + N*DMZ1
         IMPHELL(3*A,3*A-1)= IMPHELL(3*A,3*A-1) + M*DNZ1 + N*DMZ1
C
C DX1DX4,DX4DX1
C
         N=(V*DUX1 - U*DVX1)/V**2
         VIXEN3= (1.0D0+XE*DLAMBDAX1)*(1.0D0+XE*DMUX4) + (YE**2 + ZE**2)*DLAMBDAX1*DMUX4
         DNX4= (V*VIXEN3 + DVX4*DUX1 - U*DQX4*DPX1 - DVX1*DUX4)/V**2 - 2.0D0*DVX4*N/V

         IMPHELL(3*A-2,3*D-2)= IMPHELL(3*A-2,3*D-2) + M*DNX4 + N*DMX4
         IMPHELL(3*D-2,3*A-2)= IMPHELL(3*D-2,3*A-2) + M*DNX4 + N*DMX4
C
C DY1DY4,DY4DY1
C
         N=(V*DUY1 - U*DVY1)/V**2
         VIXEN3= (1.0D0+YE*DLAMBDAY1)*(1.0D0+YE*DMUY4) + (XE**2 + ZE**2)*DLAMBDAY1*DMUY4
         DNY4= (V*VIXEN3 + DVY4*DUY1 - U*DQY4*DPY1 - DVY1*DUY4)/V**2 - 2.0D0*DVY4*N/V

         IMPHELL(3*A-1,3*D-1)= IMPHELL(3*A-1,3*D-1) + M*DNY4 + N*DMY4
         IMPHELL(3*D-1,3*A-1)= IMPHELL(3*D-1,3*A-1) + M*DNY4 + N*DMY4
C
C DZ1DZ4,DZ4DZ1
C
         N=(V*DUZ1 - U*DVZ1)/V**2
         VIXEN3= (1.0D0+ZE*DLAMBDAZ1)*(1.0D0+ZE*DMUZ4) + (YE**2 + XE**2)*DLAMBDAZ1*DMUZ4
         DNZ4= (V*VIXEN3 + DVZ4*DUZ1 - U*DQZ4*DPZ1 - DVZ1*DUZ4)/V**2 - 2.0D0*DVZ4*N/V

         IMPHELL(3*A,3*D)= IMPHELL(3*A,3*D) + M*DNZ4 + N*DMZ4
         IMPHELL(3*D,3*A)= IMPHELL(3*D,3*A) + M*DNZ4 + N*DMZ4
C
C DX1DY4,DY4DX1
C
         N=(V*DUX1 - U*DVX1)/V**2
         VIXEN3= (1.0D0+XE*DLAMBDAX1)*XE*DMUY4 + YE*DLAMBDAX1*(1.0D0+YE*DMUY4) 
     1          + ZE**2*DLAMBDAX1*DMUY4
         DNY4= (V*VIXEN3 + DVY4*DUX1 - U*DQY4*DPX1 - DVX1*DUY4)/V**2 - 2.0D0*DVY4*N/V

         IMPHELL(3*A-2,3*D-1)= IMPHELL(3*A-2,3*D-1) + M*DNY4 + N*DMY4
         IMPHELL(3*D-1,3*A-2)= IMPHELL(3*D-1,3*A-2) + M*DNY4 + N*DMY4
C
C DX1DZ4,DZ4DX1
C
         N=(V*DUX1 - U*DVX1)/V**2
         VIXEN3= (1.0D0+XE*DLAMBDAX1)*XE*DMUZ4 + ZE*DLAMBDAX1*(1.0D0+ZE*DMUZ4)
     1          + YE**2*DLAMBDAX1*DMUZ4
         DNZ4= (V*VIXEN3 + DVZ4*DUX1 - U*DQZ4*DPX1 - DVX1*DUZ4)/V**2 - 2.0D0*DVZ4*N/V

         IMPHELL(3*A-2,3*D)= IMPHELL(3*A-2,3*D) + M*DNZ4 + N*DMZ4
         IMPHELL(3*D,3*A-2)= IMPHELL(3*D,3*A-2) + M*DNZ4 + N*DMZ4
C
C DY1DZ4,DZ4DY1
C
         N=(V*DUY1 - U*DVY1)/V**2
         VIXEN3= (1.0D0+YE*DLAMBDAY1)*YE*DMUZ4 + ZE*DLAMBDAY1*(1.0D0+ZE*DMUZ4)
     1          + XE**2*DLAMBDAY1*DMUZ4
         DNZ4= (V*VIXEN3 + DVZ4*DUY1 - U*DQZ4*DPY1 - DVY1*DUZ4)/V**2 - 2.0D0*DVZ4*N/V

         IMPHELL(3*A-1,3*D)= IMPHELL(3*A-1,3*D) + M*DNZ4 + N*DMZ4
         IMPHELL(3*D,3*A-1)= IMPHELL(3*D,3*A-1) + M*DNZ4 + N*DMZ4
C
C DY1DX4,DX4DY1
C
         N=(V*DUY1 - U*DVY1)/V**2
         VIXEN3= (1.0D0+YE*DLAMBDAY1)*YE*DMUX4 + XE*DLAMBDAY1*(1.0D0+XE*DMUX4)
     1          + ZE**2*DLAMBDAY1*DMUX4
         DNX4= (V*VIXEN3 + DVX4*DUY1 - U*DQX4*DPY1 - DVY1*DUX4)/V**2 - 2.0D0*DVX4*N/V

         IMPHELL(3*A-1,3*D-2)= IMPHELL(3*A-1,3*D-2) + M*DNX4 + N*DMX4
         IMPHELL(3*D-2,3*A-1)= IMPHELL(3*D-2,3*A-1) + M*DNX4 + N*DMX4
C
C DZ1DX4,DX4DZ1
C
         N=(V*DUZ1 - U*DVZ1)/V**2
         VIXEN3= (1.0D0+ZE*DLAMBDAZ1)*ZE*DMUX4 + XE*DLAMBDAZ1*(1.0D0+XE*DMUX4)
     1          + YE**2*DLAMBDAZ1*DMUX4
         DNX4= (V*VIXEN3 + DVX4*DUZ1 - U*DQX4*DPZ1 - DVZ1*DUX4)/V**2 - 2.0D0*DVX4*N/V

         IMPHELL(3*A,3*D-2)= IMPHELL(3*A,3*D-2) + M*DNX4 + N*DMX4
         IMPHELL(3*D-2,3*A)= IMPHELL(3*D-2,3*A) + M*DNX4 + N*DMX4
C
C DZ1DY4,DY4DZ1
C
         N=(V*DUZ1 - U*DVZ1)/V**2
         VIXEN3= (1.0D0+ZE*DLAMBDAZ1)*ZE*DMUY4 + YE*DLAMBDAZ1*(1.0D0+YE*DMUY4)
     1          + XE**2*DLAMBDAZ1*DMUY4
         DNY4= (V*VIXEN3 + DVY4*DUZ1 - U*DQY4*DPZ1 - DVZ1*DUY4)/V**2 - 2.0D0*DVY4*N/V

         IMPHELL(3*A,3*D-1)= IMPHELL(3*A,3*D-1) + M*DNY4 + N*DMY4
         IMPHELL(3*D-1,3*A)= IMPHELL(3*D-1,3*A) + M*DNY4 + N*DMY4
C
C DX4DX4
C
         N=(V*DUX4 - U*DVX4)/V**2
         VIXEN3= (VIXEN1/VIXEN2)*((1.0D0+XE*DMUX4)**2 + (YE**2 + ZE**2)*DMUX4**2 - DQX4**2)
         DNX4= -U*VIXEN3/V**2 - 2.0D0*DVX4*N/V

         IMPHELL(3*D-2,3*D-2)= IMPHELL(3*D-2,3*D-2) + M*DNX4 + N*DMX4
C
C DY4DY4
C
         N=(V*DUY4 - U*DVY4)/V**2
         VIXEN3= (VIXEN1/VIXEN2)*((1.0D0+YE*DMUY4)**2 + (XE**2 + ZE**2)*DMUY4**2 - DQY4**2)
         DNY4= -U*VIXEN3/V**2 - 2.0D0*DVY4*N/V

         IMPHELL(3*D-1,3*D-1)= IMPHELL(3*D-1,3*D-1) + M*DNY4 + N*DMY4
C
C DZ4DZ4
C
         N=(V*DUZ4 - U*DVZ4)/V**2
         VIXEN3= (VIXEN1/VIXEN2)*((1.0D0+ZE*DMUZ4)**2 + (YE**2 + XE**2)*DMUZ4**2 - DQZ4**2)
         DNZ4= -U*VIXEN3/V**2 - 2.0D0*DVZ4*N/V

         IMPHELL(3*D,3*D)= IMPHELL(3*D,3*D) + M*DNZ4 + N*DMZ4
C
C DX4DY4,DY4DX4
C
         N=(V*DUX4 - U*DVX4)/V**2
         VIXEN3= (VIXEN1/VIXEN2)*((1.0D0+XE*DMUX4)*XE*DMUY4 + YE*DMUX4*(1.0D0+YE*DMUY4)
     1          + ZE**2*DMUX4*DMUY4 - DQX4*DQY4)
         DNY4= (DVY4*DUX4 - U*VIXEN3 - DUY4*DVX4)/V**2 - 2.0D0*DVY4*N/V

         IMPHELL(3*D-2,3*D-1)= IMPHELL(3*D-2,3*D-1) + M*DNY4 + N*DMY4
         IMPHELL(3*D-1,3*D-2)= IMPHELL(3*D-1,3*D-2) + M*DNY4 + N*DMY4
C
C DX4DZ4,DZ4DX4
C
         N=(V*DUX4 - U*DVX4)/V**2
         VIXEN3= (VIXEN1/VIXEN2)*((1.0D0+XE*DMUX4)*XE*DMUZ4 + ZE*DMUX4*(1.0D0+ZE*DMUZ4)
     1          + YE**2*DMUX4*DMUZ4 - DQX4*DQZ4)
         DNZ4= (DVZ4*DUX4 - U*VIXEN3 - DUZ4*DVX4)/V**2 - 2.0D0*DVZ4*N/V

         IMPHELL(3*D-2,3*D)= IMPHELL(3*D-2,3*D) + M*DNZ4 + N*DMZ4
         IMPHELL(3*D,3*D-2)= IMPHELL(3*D,3*D-2) + M*DNZ4 + N*DMZ4
C
C DY4DZ4,DZ4DY4
C
         N=(V*DUY4 - U*DVY4)/V**2
         VIXEN3= (VIXEN1/VIXEN2)*((1.0D0+YE*DMUY4)*YE*DMUZ4 + ZE*DMUY4*(1.0D0+ZE*DMUZ4)
     1          + XE**2*DMUY4*DMUZ4 - DQY4*DQZ4)
         DNZ4= (DVZ4*DUY4 - U*VIXEN3 - DUZ4*DVY4)/V**2 - 2.0D0*DVZ4*N/V

         IMPHELL(3*D-1,3*D)= IMPHELL(3*D-1,3*D) + M*DNZ4 + N*DMZ4
         IMPHELL(3*D,3*D-1)= IMPHELL(3*D,3*D-1) + M*DNZ4 + N*DMZ4
C
C DX1DX2,DX2DX1
C
         N=(V*DUX1 - U*DVX1)/V**2
         VIXEN4= 1.0D0/RBC2 - 2.0D0*(X(B)-X(C))**2/RBC4
         VIXEN5= XP*(VIXEN4*XE-DLAMBDAX1) + (1.0D0+XE*DLAMBDAX1)*(XE*DLAMBDAX2-1.0D0-LAMBDA)
     1           + YP*YE*VIXEN4 + YE**2*DLAMBDAX1*DLAMBDAX2 + ZP*ZE*VIXEN4 
     2           + ZE**2*DLAMBDAX1*DLAMBDAX2 - DPX1*DPX2
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN5 + DPX1*DQX2 
         VIXEN6= XPP*(VIXEN4*XE-DLAMBDAX1) + (1.0D0+XE*DLAMBDAX1)*(XE*DMUX2-1.0D0-MU) +YPP*YE*VIXEN4
     1           + YE**2*DLAMBDAX1*DMUX2 + ZPP*ZE*VIXEN4 + ZE**2*DLAMBDAX1*DMUX2
         DNX2= (V*VIXEN6 + DVX2*DUX1 - U*VIXEN3 - DUX2*DVX1)/V**2 - 2.0D0*DVX2*N/V

         IMPHELL(3*A-2,3*B-2)= IMPHELL(3*A-2,3*B-2)+ M*DNX2 + N*DMX2
         IMPHELL(3*B-2,3*A-2)= IMPHELL(3*B-2,3*A-2)+ M*DNX2 + N*DMX2
C
C DY1DY2,DY2DY1
C
         N=(V*DUY1 - U*DVY1)/V**2
         VIXEN4= 1.0D0/RBC2 - 2.0D0*(Y(B)-Y(C))**2/RBC4
         VIXEN5= YP*(VIXEN4*YE-DLAMBDAY1) + (1.0D0+YE*DLAMBDAY1)*(YE*DLAMBDAY2-1.0D0-LAMBDA)
     1           + XP*XE*VIXEN4 + XE**2*DLAMBDAY1*DLAMBDAY2 + ZP*ZE*VIXEN4
     2           + ZE**2*DLAMBDAY1*DLAMBDAY2 - DPY1*DPY2
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN5 + DPY1*DQY2
         VIXEN6= YPP*(VIXEN4*YE-DLAMBDAY1) + (1.0D0+YE*DLAMBDAY1)*(YE*DMUY2-1.0D0-MU) +XPP*XE*VIXEN4
     1           + XE**2*DLAMBDAY1*DMUY2 + ZPP*ZE*VIXEN4 + ZE**2*DLAMBDAY1*DMUY2
         DNY2= (V*VIXEN6 + DVY2*DUY1 - U*VIXEN3 - DUY2*DVY1)/V**2 - 2.0D0*DVY2*N/V

         IMPHELL(3*A-1,3*B-1)= IMPHELL(3*A-1,3*B-1)+ M*DNY2 + N*DMY2
         IMPHELL(3*B-1,3*A-1)= IMPHELL(3*B-1,3*A-1)+ M*DNY2 + N*DMY2
C
C DZ1DZ2,DZ2DZ1
C
         N=(V*DUZ1 - U*DVZ1)/V**2
         VIXEN4= 1.0D0/RBC2 - 2.0D0*(Z(B)-Z(C))**2/RBC4
         VIXEN5= ZP*(VIXEN4*ZE-DLAMBDAZ1) + (1.0D0+ZE*DLAMBDAZ1)*(ZE*DLAMBDAZ2-1.0D0-LAMBDA)
     1           + XP*XE*VIXEN4 + XE**2*DLAMBDAZ1*DLAMBDAZ2 + YP*YE*VIXEN4
     2           + YE**2*DLAMBDAZ1*DLAMBDAZ2 - DPZ1*DPZ2
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN5 + DPZ1*DQZ2
         VIXEN6= ZPP*(VIXEN4*ZE-DLAMBDAZ1) + (1.0D0+ZE*DLAMBDAZ1)*(ZE*DMUZ2-1.0D0-MU) +XPP*XE*VIXEN4
     1           + XE**2*DLAMBDAZ1*DMUZ2 + YPP*YE*VIXEN4 + YE**2*DLAMBDAZ1*DMUZ2
         DNZ2= (V*VIXEN6 + DVZ2*DUZ1 - U*VIXEN3 - DUZ2*DVZ1)/V**2 - 2.0D0*DVZ2*N/V

         IMPHELL(3*A,3*B)= IMPHELL(3*A,3*B)+ M*DNZ2 + N*DMZ2
         IMPHELL(3*B,3*A)= IMPHELL(3*B,3*A)+ M*DNZ2 + N*DMZ2
C
C DX1DY2,DY2DX1
C
         N=(V*DUX1 - U*DVX1)/V**2
         VIXEN4= -2.0D0*(X(B)-X(C))*(Y(B)-Y(C))/RBC4
         VIXEN5= XP*XE*VIXEN4 + (1.0D0+XE*DLAMBDAX1)*XE*DLAMBDAY2 + YP*(YE*VIXEN4-DLAMBDAX1)
     1          + YE*DLAMBDAX1*(YE*DLAMBDAY2-1.0D0-LAMBDA) + ZP*ZE*VIXEN4 
     2          + ZE**2*DLAMBDAX1*DLAMBDAY2 - DPX1*DPY2
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN5 + DPX1*DQY2
         VIXEN6= XPP*XE*VIXEN4 + (1.0D0+XE*DLAMBDAX1)*XE*DMUY2 + YPP*(YE*VIXEN4-DLAMBDAX1)
     1          + YE*DLAMBDAX1*(YE*DMUY2-1.0D0-MU) + ZPP*ZE*VIXEN4 + ZE**2*DLAMBDAX1*DMUY2
         DNY2= (V*VIXEN6 + DVY2*DUX1 - U*VIXEN3 - DUY2*DVX1)/V**2 - 2.0D0*DVY2*N/V

         IMPHELL(3*A-2,3*B-1)= IMPHELL(3*A-2,3*B-1)+ M*DNY2 + N*DMY2
         IMPHELL(3*B-1,3*A-2)= IMPHELL(3*B-1,3*A-2)+ M*DNY2 + N*DMY2
C
C DX1DZ2,DZ2DX1
C
         N=(V*DUX1 - U*DVX1)/V**2
         VIXEN4= -2.0D0*(X(B)-X(C))*(Z(B)-Z(C))/RBC4
         VIXEN5= XP*XE*VIXEN4 + (1.0D0+XE*DLAMBDAX1)*XE*DLAMBDAZ2 + ZP*(ZE*VIXEN4-DLAMBDAX1)
     1          + ZE*DLAMBDAX1*(ZE*DLAMBDAZ2-1.0D0-LAMBDA) + YP*YE*VIXEN4
     2          + YE**2*DLAMBDAX1*DLAMBDAZ2 - DPX1*DPZ2
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN5 + DPX1*DQZ2
         VIXEN6= XPP*XE*VIXEN4 + (1.0D0+XE*DLAMBDAX1)*XE*DMUZ2 + ZPP*(ZE*VIXEN4-DLAMBDAX1)
     1          + ZE*DLAMBDAX1*(ZE*DMUZ2-1.0D0-MU) + YPP*YE*VIXEN4 + YE**2*DLAMBDAX1*DMUZ2
         DNZ2= (V*VIXEN6 + DVZ2*DUX1 - U*VIXEN3 - DUZ2*DVX1)/V**2 - 2.0D0*DVZ2*N/V

         IMPHELL(3*A-2,3*B)= IMPHELL(3*A-2,3*B)+ M*DNZ2 + N*DMZ2
         IMPHELL(3*B,3*A-2)= IMPHELL(3*B,3*A-2)+ M*DNZ2 + N*DMZ2
C
C DY1DZ2,DZ2DY1
C
         N=(V*DUY1 - U*DVY1)/V**2
         VIXEN4= -2.0D0*(Y(B)-Y(C))*(Z(B)-Z(C))/RBC4
         VIXEN5= YP*YE*VIXEN4 + (1.0D0+YE*DLAMBDAY1)*YE*DLAMBDAZ2 + ZP*(ZE*VIXEN4-DLAMBDAY1)
     1          + ZE*DLAMBDAY1*(ZE*DLAMBDAZ2-1.0D0-LAMBDA) + XP*XE*VIXEN4
     2          + XE**2*DLAMBDAY1*DLAMBDAZ2 - DPY1*DPZ2
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN5 + DPY1*DQZ2
         VIXEN6= YPP*YE*VIXEN4 + (1.0D0+YE*DLAMBDAY1)*YE*DMUZ2 + ZPP*(ZE*VIXEN4-DLAMBDAY1)
     1          + ZE*DLAMBDAY1*(ZE*DMUZ2-1.0D0-MU) + XPP*XE*VIXEN4 + XE**2*DLAMBDAY1*DMUZ2
         DNZ2= (V*VIXEN6 + DVZ2*DUY1 - U*VIXEN3 - DUZ2*DVY1)/V**2 - 2.0D0*DVZ2*N/V

         IMPHELL(3*A-1,3*B)= IMPHELL(3*A-1,3*B)+ M*DNZ2 + N*DMZ2
         IMPHELL(3*B,3*A-1)= IMPHELL(3*B,3*A-1)+ M*DNZ2 + N*DMZ2
C
C DY1DX2,DX2DY1
C
         N=(V*DUY1 - U*DVY1)/V**2
         VIXEN4= -2.0D0*(Y(B)-Y(C))*(X(B)-X(C))/RBC4
         VIXEN5= YP*YE*VIXEN4 + (1.0D0+YE*DLAMBDAY1)*YE*DLAMBDAX2 + XP*(XE*VIXEN4-DLAMBDAY1)
     1          + XE*DLAMBDAY1*(XE*DLAMBDAX2-1.0D0-LAMBDA) + ZP*ZE*VIXEN4
     2          + ZE**2*DLAMBDAY1*DLAMBDAX2 - DPY1*DPX2
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN5 + DPY1*DQX2
         VIXEN6= YPP*YE*VIXEN4 + (1.0D0+YE*DLAMBDAY1)*YE*DMUX2 + XPP*(XE*VIXEN4-DLAMBDAY1)
     1          + XE*DLAMBDAY1*(XE*DMUX2-1.0D0-MU) + ZPP*ZE*VIXEN4 + ZE**2*DLAMBDAY1*DMUX2
         DNX2= (V*VIXEN6 + DVX2*DUY1 - U*VIXEN3 - DUX2*DVY1)/V**2 - 2.0D0*DVX2*N/V

         IMPHELL(3*A-1,3*B-2)= IMPHELL(3*A-1,3*B-2)+ M*DNX2 + N*DMX2
         IMPHELL(3*B-2,3*A-1)= IMPHELL(3*B-2,3*A-1)+ M*DNX2 + N*DMX2
C
C DZ1DX2,DX2DZ1
C
         N=(V*DUZ1 - U*DVZ1)/V**2
         VIXEN4= -2.0D0*(Z(B)-Z(C))*(X(B)-X(C))/RBC4
         VIXEN5= ZP*ZE*VIXEN4 + (1.0D0+ZE*DLAMBDAZ1)*ZE*DLAMBDAX2 + XP*(XE*VIXEN4-DLAMBDAZ1)
     1          + XE*DLAMBDAZ1*(XE*DLAMBDAX2-1.0D0-LAMBDA) + YP*YE*VIXEN4
     2          + YE**2*DLAMBDAZ1*DLAMBDAX2 - DPZ1*DPX2
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN5 + DPZ1*DQX2
         VIXEN6= ZPP*ZE*VIXEN4 + (1.0D0+ZE*DLAMBDAZ1)*ZE*DMUX2 + XPP*(XE*VIXEN4-DLAMBDAZ1)
     1          + XE*DLAMBDAZ1*(XE*DMUX2-1.0D0-MU) + YPP*YE*VIXEN4 + YE**2*DLAMBDAZ1*DMUX2
         DNX2= (V*VIXEN6 + DVX2*DUZ1 - U*VIXEN3 - DUX2*DVZ1)/V**2 - 2.0D0*DVX2*N/V

         IMPHELL(3*A,3*B-2)= IMPHELL(3*A,3*B-2)+ M*DNX2 + N*DMX2
         IMPHELL(3*B-2,3*A)= IMPHELL(3*B-2,3*A)+ M*DNX2 + N*DMX2
C
C DZ1DY2,DY2DZ1
C
         N=(V*DUZ1 - U*DVZ1)/V**2
         VIXEN4= -2.0D0*(Z(B)-Z(C))*(Y(B)-Y(C))/RBC4
         VIXEN5= ZP*ZE*VIXEN4 + (1.0D0+ZE*DLAMBDAZ1)*ZE*DLAMBDAY2 + YP*(YE*VIXEN4-DLAMBDAZ1)
     1          + YE*DLAMBDAZ1*(YE*DLAMBDAY2-1.0D0-LAMBDA) + XP*XE*VIXEN4
     2          + XE**2*DLAMBDAZ1*DLAMBDAY2 - DPZ1*DPY2
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN5 + DPZ1*DQY2
         VIXEN6= ZPP*ZE*VIXEN4 + (1.0D0+ZE*DLAMBDAZ1)*ZE*DMUY2 + YPP*(YE*VIXEN4-DLAMBDAZ1)
     1          + YE*DLAMBDAZ1*(YE*DMUY2-1.0D0-MU) + XPP*XE*VIXEN4 + XE**2*DLAMBDAZ1*DMUY2
         DNY2= (V*VIXEN6 + DVY2*DUZ1 - U*VIXEN3 - DUY2*DVZ1)/V**2 - 2.0D0*DVY2*N/V

         IMPHELL(3*A,3*B-1)= IMPHELL(3*A,3*B-1)+ M*DNY2 + N*DMY2
         IMPHELL(3*B-1,3*A)= IMPHELL(3*B-1,3*A)+ M*DNY2 + N*DMY2
C
C DX1DX3,DX3DX1
C
         N=(V*DUX1 - U*DVX1)/V**2
         VIXEN4= 2.0D0*(X(B)-X(C))**2/RBC4 - 1.0D0/RBC2
         VIXEN5= XP*(XE*VIXEN4+DLAMBDAX1) + (1.0D0+XE*DLAMBDAX1)*(XE*DLAMBDAX3+LAMBDA)
     1          + (YP*YE+ZP*ZE)*VIXEN4 + (YE**2+ZE**2)*DLAMBDAX1*DLAMBDAX3 - DPX1*DPX3
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN5 + DPX1*DQX3
         VIXEN6= XPP*(XE*VIXEN4+DLAMBDAX1) + (1.0D0+XE*DLAMBDAX1)*(XE*DMUX3+MU)
     1          + (YPP*YE+ZPP*ZE)*VIXEN4 + (YE**2+ZE**2)*DLAMBDAX1*DMUX3
         DNX3= (V*VIXEN6 + DVX3*DUX1 - U*VIXEN3 - DUX3*DVX1)/V**2 - 2.0D0*DVX3*N/V

         IMPHELL(3*A-2,3*C-2)= IMPHELL(3*A-2,3*C-2)+ M*DNX3 + N*DMX3
         IMPHELL(3*C-2,3*A-2)= IMPHELL(3*C-2,3*A-2)+ M*DNX3 + N*DMX3
C
C DY1DY3,DY3DY1
C
         N=(V*DUY1 - U*DVY1)/V**2
         VIXEN4= 2.0D0*(Y(B)-Y(C))**2/RBC4 - 1.0D0/RBC2
         VIXEN5= YP*(YE*VIXEN4+DLAMBDAY1) + (1.0D0+YE*DLAMBDAY1)*(YE*DLAMBDAY3+LAMBDA)
     1          + (XP*XE+ZP*ZE)*VIXEN4 + (XE**2+ZE**2)*DLAMBDAY1*DLAMBDAY3 - DPY1*DPY3
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN5 + DPY1*DQY3
         VIXEN6= YPP*(YE*VIXEN4+DLAMBDAY1) + (1.0D0+YE*DLAMBDAY1)*(YE*DMUY3+MU)
     1          + (XPP*XE+ZPP*ZE)*VIXEN4 + (XE**2+ZE**2)*DLAMBDAY1*DMUY3
         DNY3= (V*VIXEN6 + DVY3*DUY1 - U*VIXEN3 - DUY3*DVY1)/V**2 - 2.0D0*DVY3*N/V

         IMPHELL(3*A-1,3*C-1)= IMPHELL(3*A-1,3*C-1)+ M*DNY3 + N*DMY3
         IMPHELL(3*C-1,3*A-1)= IMPHELL(3*C-1,3*A-1)+ M*DNY3 + N*DMY3
C
C DZ1DZ3,DZ3DZ1
C
         N=(V*DUZ1 - U*DVZ1)/V**2
         VIXEN4= 2.0D0*(Z(B)-Z(C))**2/RBC4 - 1.0D0/RBC2
         VIXEN5= ZP*(ZE*VIXEN4+DLAMBDAZ1) + (1.0D0+ZE*DLAMBDAZ1)*(ZE*DLAMBDAZ3+LAMBDA)
     1          + (XP*XE+YP*YE)*VIXEN4 + (XE**2+YE**2)*DLAMBDAZ1*DLAMBDAZ3 - DPZ1*DPZ3
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN5 + DPZ1*DQZ3
         VIXEN6= ZPP*(ZE*VIXEN4+DLAMBDAZ1) + (1.0D0+ZE*DLAMBDAZ1)*(ZE*DMUZ3+MU)
     1          + (XPP*XE+YPP*YE)*VIXEN4 + (XE**2+YE**2)*DLAMBDAZ1*DMUZ3
         DNZ3= (V*VIXEN6 + DVZ3*DUZ1 - U*VIXEN3 - DUZ3*DVZ1)/V**2 - 2.0D0*DVZ3*N/V

         IMPHELL(3*A,3*C)= IMPHELL(3*A,3*C)+ M*DNZ3 + N*DMZ3
         IMPHELL(3*C,3*A)= IMPHELL(3*C,3*A)+ M*DNZ3 + N*DMZ3
C
C DX1DY3
C
         N=(V*DUX1 - U*DVX1)/V**2
         VIXEN4= 2.0D0*(X(B)-X(C))*(Y(B)-Y(C))/RBC4
         VIXEN5= XP*XE*VIXEN4 + (1.0D0+XE*DLAMBDAX1)*XE*DLAMBDAY3 
     1          + YP*(YE*VIXEN4+DLAMBDAX1) + YE*DLAMBDAX1*(YE*DLAMBDAY3+LAMBDA) 
     2          + ZP*ZE*VIXEN4 + ZE**2*DLAMBDAX1*DLAMBDAY3 - DPX1*DPY3
         VIXEN3=(VIXEN2/VIXEN1)*VIXEN5 + DPX1*DQY3
         VIXEN6= XPP*XE*VIXEN4 + (1.0D0+XE*DLAMBDAX1)*XE*DMUY3 + YPP*(YE*VIXEN4+DLAMBDAX1)
     1          + YE*DLAMBDAX1*(YE*DMUY3+MU) + ZPP*ZE*VIXEN4 + ZE**2*DLAMBDAX1*DMUY3
         DNY3= (V*VIXEN6 + DVY3*DUX1 - U*VIXEN3 - DUY3*DVX1)/V**2 - 2.0D0*DVY3*N/V

         IMPHELL(3*A-2,3*C-1)= IMPHELL(3*A-2,3*C-1)+ M*DNY3 + N*DMY3
         IMPHELL(3*C-1,3*A-2)= IMPHELL(3*C-1,3*A-2)+ M*DNY3 + N*DMY3
C
C DX1DZ3
C
         N=(V*DUX1 - U*DVX1)/V**2
         VIXEN4= 2.0D0*(X(B)-X(C))*(Z(B)-Z(C))/RBC4
         VIXEN5= XP*XE*VIXEN4 + (1.0D0+XE*DLAMBDAX1)*XE*DLAMBDAZ3
     1          + ZP*(ZE*VIXEN4+DLAMBDAX1) + ZE*DLAMBDAX1*(ZE*DLAMBDAZ3+LAMBDA)
     2          + YP*YE*VIXEN4 + YE**2*DLAMBDAX1*DLAMBDAZ3 - DPX1*DPZ3
         VIXEN3=(VIXEN2/VIXEN1)*VIXEN5 + DPX1*DQZ3
         VIXEN6= XPP*XE*VIXEN4 + (1.0D0+XE*DLAMBDAX1)*XE*DMUZ3 + ZPP*(ZE*VIXEN4+DLAMBDAX1)
     1          + ZE*DLAMBDAX1*(ZE*DMUZ3+MU) + YPP*YE*VIXEN4 + YE**2*DLAMBDAX1*DMUZ3
         DNZ3= (V*VIXEN6 + DVZ3*DUX1 - U*VIXEN3 - DUZ3*DVX1)/V**2 - 2.0D0*DVZ3*N/V

         IMPHELL(3*A-2,3*C)= IMPHELL(3*A-2,3*C)+ M*DNZ3 + N*DMZ3
         IMPHELL(3*C,3*A-2)= IMPHELL(3*C,3*A-2)+ M*DNZ3 + N*DMZ3
C
C DY1DZ3
C
         N=(V*DUY1 - U*DVY1)/V**2
         VIXEN4= 2.0D0*(Y(B)-Y(C))*(Z(B)-Z(C))/RBC4
         VIXEN5= YP*YE*VIXEN4 + (1.0D0+YE*DLAMBDAY1)*YE*DLAMBDAZ3
     1          + ZP*(ZE*VIXEN4+DLAMBDAY1) + ZE*DLAMBDAY1*(ZE*DLAMBDAZ3+LAMBDA)
     2          + XP*XE*VIXEN4 + XE**2*DLAMBDAY1*DLAMBDAZ3 - DPY1*DPZ3
         VIXEN3=(VIXEN2/VIXEN1)*VIXEN5 + DPY1*DQZ3
         VIXEN6= YPP*YE*VIXEN4 + (1.0D0+YE*DLAMBDAY1)*YE*DMUZ3 + ZPP*(ZE*VIXEN4+DLAMBDAY1)
     1          + ZE*DLAMBDAY1*(ZE*DMUZ3+MU) + XPP*XE*VIXEN4 + XE**2*DLAMBDAY1*DMUZ3
         DNZ3= (V*VIXEN6 + DVZ3*DUY1 - U*VIXEN3 - DUZ3*DVY1)/V**2 - 2.0D0*DVZ3*N/V

         IMPHELL(3*A-1,3*C)= IMPHELL(3*A-1,3*C)+ M*DNZ3 + N*DMZ3
         IMPHELL(3*C,3*A-1)= IMPHELL(3*C,3*A-1)+ M*DNZ3 + N*DMZ3
C
C DY1DX3
C
         N=(V*DUY1 - U*DVY1)/V**2
         VIXEN4= 2.0D0*(Y(B)-Y(C))*(X(B)-X(C))/RBC4
         VIXEN5= YP*YE*VIXEN4 + (1.0D0+YE*DLAMBDAY1)*YE*DLAMBDAX3
     1          + XP*(XE*VIXEN4+DLAMBDAY1) + XE*DLAMBDAY1*(XE*DLAMBDAX3+LAMBDA)
     2          + ZP*ZE*VIXEN4 + ZE**2*DLAMBDAY1*DLAMBDAX3 - DPY1*DPX3
         VIXEN3=(VIXEN2/VIXEN1)*VIXEN5 + DPY1*DQX3
         VIXEN6= YPP*YE*VIXEN4 + (1.0D0+YE*DLAMBDAY1)*YE*DMUX3 + XPP*(XE*VIXEN4+DLAMBDAY1)
     1          + XE*DLAMBDAY1*(XE*DMUX3+MU) + ZPP*ZE*VIXEN4 + ZE**2*DLAMBDAY1*DMUX3
         DNX3= (V*VIXEN6 + DVX3*DUY1 - U*VIXEN3 - DUX3*DVY1)/V**2 - 2.0D0*DVX3*N/V

         IMPHELL(3*A-1,3*C-2)= IMPHELL(3*A-1,3*C-2)+ M*DNX3 + N*DMX3
         IMPHELL(3*C-2,3*A-1)= IMPHELL(3*C-2,3*A-1)+ M*DNX3 + N*DMX3
C
C DZ1DX3
C
         N=(V*DUZ1 - U*DVZ1)/V**2
         VIXEN4= 2.0D0*(Z(B)-Z(C))*(X(B)-X(C))/RBC4
         VIXEN5= ZP*ZE*VIXEN4 + (1.0D0+ZE*DLAMBDAZ1)*ZE*DLAMBDAX3
     1          + XP*(XE*VIXEN4+DLAMBDAZ1) + XE*DLAMBDAZ1*(XE*DLAMBDAX3+LAMBDA)
     2          + YP*YE*VIXEN4 + YE**2*DLAMBDAZ1*DLAMBDAX3 - DPZ1*DPX3
         VIXEN3=(VIXEN2/VIXEN1)*VIXEN5 + DPZ1*DQX3
         VIXEN6= ZPP*ZE*VIXEN4 + (1.0D0+ZE*DLAMBDAZ1)*ZE*DMUX3 + XPP*(XE*VIXEN4+DLAMBDAZ1)
     1          + XE*DLAMBDAZ1*(XE*DMUX3+MU) + YPP*YE*VIXEN4 + YE**2*DLAMBDAZ1*DMUX3
         DNX3= (V*VIXEN6 + DVX3*DUZ1 - U*VIXEN3 - DUX3*DVZ1)/V**2 - 2.0D0*DVX3*N/V

         IMPHELL(3*A,3*C-2)= IMPHELL(3*A,3*C-2)+ M*DNX3 + N*DMX3
         IMPHELL(3*C-2,3*A)= IMPHELL(3*C-2,3*A)+ M*DNX3 + N*DMX3
C
C DZ1DY3
C
         N=(V*DUZ1 - U*DVZ1)/V**2
         VIXEN4= 2.0D0*(Z(B)-Z(C))*(Y(B)-Y(C))/RBC4
         VIXEN5= ZP*ZE*VIXEN4 + (1.0D0+ZE*DLAMBDAZ1)*ZE*DLAMBDAY3
     1          + YP*(YE*VIXEN4+DLAMBDAZ1) + YE*DLAMBDAZ1*(YE*DLAMBDAY3+LAMBDA)
     2          + XP*XE*VIXEN4 + XE**2*DLAMBDAZ1*DLAMBDAY3 - DPZ1*DPY3
         VIXEN3=(VIXEN2/VIXEN1)*VIXEN5 + DPZ1*DQY3
         VIXEN6= ZPP*ZE*VIXEN4 + (1.0D0+ZE*DLAMBDAZ1)*ZE*DMUY3 + YPP*(YE*VIXEN4+DLAMBDAZ1)
     1          + YE*DLAMBDAZ1*(YE*DMUY3+MU) + XPP*XE*VIXEN4 + XE**2*DLAMBDAZ1*DMUY3
         DNY3= (V*VIXEN6 + DVY3*DUZ1 - U*VIXEN3 - DUY3*DVZ1)/V**2 - 2.0D0*DVY3*N/V

         IMPHELL(3*A,3*C-1)= IMPHELL(3*A,3*C-1)+ M*DNY3 + N*DMY3
         IMPHELL(3*C-1,3*A)= IMPHELL(3*C-1,3*A)+ M*DNY3 + N*DMY3
C
C DX4DX2,DX2DX4
C
         N=(V*DUX4 - U*DVX4)/V**2
         VIXEN4= 1.0D0/RBC2 - 2.0D0*(X(B)-X(C))**2/RBC4
         VIXEN5= XPP*(XE*VIXEN4-DMUX4) + (1.0D0+XE*DMUX4)
     1          *(XE*DMUX2-1.0D0-MU) + YPP*YE*VIXEN4
     2          + YE**2*DMUX4*DMUX2 + ZPP*ZE*VIXEN4
     3          + ZE**2*DMUX4*DMUX2 - DQX4*DQX2
         VIXEN3=(VIXEN1/VIXEN2)*VIXEN5 + DPX2*DQX4
         VIXEN6= XP*(XE*VIXEN4-DMUX4) + (1.0D0+XE*DMUX4)
     1          *(XE*DLAMBDAX2-1.0D0-LAMBDA) + YP*YE*VIXEN4
     2          + YE**2*DMUX4*DLAMBDAX2 + ZP*ZE*VIXEN4
     3          + ZE**2*DMUX4*DLAMBDAX2
         DNX2= (V*VIXEN6 + DVX2*DUX4 - U*VIXEN3 - DUX2*DVX4)/V**2
     1        - 2.0D0*DVX2*(V*DUX4 - U*DVX4)/V**3

         IMPHELL(3*D-2,3*B-2)= IMPHELL(3*D-2,3*B-2)+ M*DNX2 + N*DMX2
         IMPHELL(3*B-2,3*D-2)= IMPHELL(3*B-2,3*D-2)+ M*DNX2 + N*DMX2
C
C DY4DY2,DY2DY4
C
         N=(V*DUY4 - U*DVY4)/V**2
         VIXEN4= 1.0D0/RBC2 - 2.0D0*(Y(B)-Y(C))**2/RBC4
         VIXEN5= YPP*(YE*VIXEN4-DMUY4) + (1.0D0+YE*DMUY4)
     1          *(YE*DMUY2-1.0D0-MU) + XPP*XE*VIXEN4
     2          + XE**2*DMUY4*DMUY2 + ZPP*ZE*VIXEN4
     3          + ZE**2*DMUY4*DMUY2 - DQY4*DQY2
         VIXEN3=(VIXEN1/VIXEN2)*VIXEN5 + DPY2*DQY4
         VIXEN6= YP*(YE*VIXEN4-DMUY4) + (1.0D0+YE*DMUY4)
     1          *(YE*DLAMBDAY2-1.0D0-LAMBDA) + XP*XE*VIXEN4
     2          + XE**2*DMUY4*DLAMBDAY2 + ZP*ZE*VIXEN4
     3          + ZE**2*DMUY4*DLAMBDAY2
         DNY2= (V*VIXEN6 + DVY2*DUY4 - U*VIXEN3 - DUY2*DVY4)/V**2
     1        - 2.0D0*DVY2*(V*DUY4 - U*DVY4)/V**3

         IMPHELL(3*D-1,3*B-1)= IMPHELL(3*D-1,3*B-1)+ M*DNY2 + N*DMY2
         IMPHELL(3*B-1,3*D-1)= IMPHELL(3*B-1,3*D-1)+ M*DNY2 + N*DMY2
C
C DZ4DZ2,DZ2DZ4
C
         N=(V*DUZ4 - U*DVZ4)/V**2
         VIXEN4= 1.0D0/RBC2 - 2.0D0*(Z(B)-Z(C))**2/RBC4
         VIXEN5= ZPP*(ZE*VIXEN4-DMUZ4) + (1.0D0+ZE*DMUZ4)
     1          *(ZE*DMUZ2-1.0D0-MU) + XPP*XE*VIXEN4
     2          + XE**2*DMUZ4*DMUZ2 + YPP*YE*VIXEN4
     3          + YE**2*DMUZ4*DMUZ2 - DQZ4*DQZ2
         VIXEN3=(VIXEN1/VIXEN2)*VIXEN5 + DPZ2*DQZ4
         VIXEN6= ZP*(ZE*VIXEN4-DMUZ4) + (1.0D0+ZE*DMUZ4)
     1          *(ZE*DLAMBDAZ2-1.0D0-LAMBDA) + XP*XE*VIXEN4
     2          + XE**2*DMUZ4*DLAMBDAZ2 + YP*YE*VIXEN4
     3          + YE**2*DMUZ4*DLAMBDAZ2
         DNZ2= (V*VIXEN6 + DVZ2*DUZ4 - U*VIXEN3 - DUZ2*DVZ4)/V**2
     1        - 2.0D0*DVZ2*(V*DUZ4 - U*DVZ4)/V**3

         IMPHELL(3*D,3*B)= IMPHELL(3*D,3*B)+ M*DNZ2 + N*DMZ2
         IMPHELL(3*B,3*D)= IMPHELL(3*B,3*D)+ M*DNZ2 + N*DMZ2
C
C DX4DY2,DY2DX4
C
         N=(V*DUX4 - U*DVX4)/V**2
         VIXEN4= -2.0D0*(X(B)-X(C))*(Y(B)-Y(C))/RBC4
         VIXEN5= XPP*XE*VIXEN4 + (1.0D0+XE*DMUX4)*XE*DMUY2 + YPP*(YE*VIXEN4-DMUX4)
     1          + YE*DMUX4*(YE*DMUY2-1.0D0-MU) + ZPP*ZE*VIXEN4 + ZE**2*DMUX4*DMUY2 - DQX4*DQY2
         VIXEN3=(VIXEN1/VIXEN2)*VIXEN5 + DPY2*DQX4
         VIXEN6= XP*XE*VIXEN4 + (1.0D0+XE*DMUX4)*XE*DLAMBDAY2 + YP*(YE*VIXEN4-DMUX4)
     1          + YE*DMUX4*(YE*DLAMBDAY2-1.0D0-LAMBDA) + ZP*ZE*VIXEN4 + ZE**2*DMUX4*DLAMBDAY2
         DNY2= (V*VIXEN6 + DVY2*DUX4 - U*VIXEN3 - DUY2*DVX4)/V**2 - 2.0D0*DVY2*N/V

         IMPHELL(3*D-2,3*B-1)= IMPHELL(3*D-2,3*B-1)+ M*DNY2 + N*DMY2
         IMPHELL(3*B-1,3*D-2)= IMPHELL(3*B-1,3*D-2)+ M*DNY2 + N*DMY2
C
C DX4DZ2,DZ2DX4
C
         N=(V*DUX4 - U*DVX4)/V**2
         VIXEN4= -2.0D0*(X(B)-X(C))*(Z(B)-Z(C))/RBC4
         VIXEN5= XPP*XE*VIXEN4 + (1.0D0+XE*DMUX4)*XE*DMUZ2 + ZPP*(ZE*VIXEN4-DMUX4)
     1          + ZE*DMUX4*(ZE*DMUZ2-1.0D0-MU) + YPP*YE*VIXEN4 + YE**2*DMUX4*DMUZ2 - DQX4*DQZ2
         VIXEN3=(VIXEN1/VIXEN2)*VIXEN5 + DPZ2*DQX4
         VIXEN6= XP*XE*VIXEN4 + (1.0D0+XE*DMUX4)*XE*DLAMBDAZ2 + ZP*(ZE*VIXEN4-DMUX4)
     1          + ZE*DMUX4*(ZE*DLAMBDAZ2-1.0D0-LAMBDA) + YP*YE*VIXEN4 + YE**2*DMUX4*DLAMBDAZ2
         DNZ2= (V*VIXEN6 + DVZ2*DUX4 - U*VIXEN3 - DUZ2*DVX4)/V**2 - 2.0D0*DVZ2*N/V

         IMPHELL(3*D-2,3*B)= IMPHELL(3*D-2,3*B)+ M*DNZ2 + N*DMZ2
         IMPHELL(3*B,3*D-2)= IMPHELL(3*B,3*D-2)+ M*DNZ2 + N*DMZ2
C
C DY4DZ2,DZ2DY4
C
         N=(V*DUY4 - U*DVY4)/V**2
         VIXEN4= -2.0D0*(Y(B)-Y(C))*(Z(B)-Z(C))/RBC4
         VIXEN5= YPP*YE*VIXEN4 + (1.0D0+YE*DMUY4)*YE*DMUZ2 + ZPP*(ZE*VIXEN4-DMUY4)
     1          + ZE*DMUY4*(ZE*DMUZ2-1.0D0-MU) + XPP*XE*VIXEN4 + XE**2*DMUY4*DMUZ2 - DQY4*DQZ2
         VIXEN3=(VIXEN1/VIXEN2)*VIXEN5 + DPZ2*DQY4
         VIXEN6= YP*YE*VIXEN4 + (1.0D0+YE*DMUY4)*YE*DLAMBDAZ2 + ZP*(ZE*VIXEN4-DMUY4)
     1          + ZE*DMUY4*(ZE*DLAMBDAZ2-1.0D0-LAMBDA) + XP*XE*VIXEN4 + XE**2*DMUY4*DLAMBDAZ2
         DNZ2= (V*VIXEN6 + DVZ2*DUY4 - U*VIXEN3 - DUZ2*DVY4)/V**2 - 2.0D0*DVZ2*N/V

         IMPHELL(3*D-1,3*B)= IMPHELL(3*D-1,3*B)+ M*DNZ2 + N*DMZ2
         IMPHELL(3*B,3*D-1)= IMPHELL(3*B,3*D-1)+ M*DNZ2 + N*DMZ2
C
C DY4DX2,DX2DY4
C
         N=(V*DUY4 - U*DVY4)/V**2
         VIXEN4= -2.0D0*(Y(B)-Y(C))*(X(B)-X(C))/RBC4
         VIXEN5= YPP*YE*VIXEN4 + (1.0D0+YE*DMUY4)*YE*DMUX2 + XPP*(XE*VIXEN4-DMUY4)
     1          + XE*DMUY4*(XE*DMUX2-1.0D0-MU) + ZPP*ZE*VIXEN4 + ZE**2*DMUY4*DMUX2 - DQY4*DQX2
         VIXEN3=(VIXEN1/VIXEN2)*VIXEN5 + DPX2*DQY4
         VIXEN6= YP*YE*VIXEN4 + (1.0D0+YE*DMUY4)*YE*DLAMBDAX2 + XP*(XE*VIXEN4-DMUY4)
     1          + XE*DMUY4*(XE*DLAMBDAX2-1.0D0-LAMBDA) + ZP*ZE*VIXEN4 + ZE**2*DMUY4*DLAMBDAX2
         DNX2= (V*VIXEN6 + DVX2*DUY4 - U*VIXEN3 - DUX2*DVY4)/V**2 - 2.0D0*DVX2*N/V

         IMPHELL(3*D-1,3*B-2)= IMPHELL(3*D-1,3*B-2)+ M*DNX2 + N*DMX2
         IMPHELL(3*B-2,3*D-1)= IMPHELL(3*B-2,3*D-1)+ M*DNX2 + N*DMX2
C
C DZ4DX2,DX2DZ4
C
         N=(V*DUZ4 - U*DVZ4)/V**2
         VIXEN4= -2.0D0*(Z(B)-Z(C))*(X(B)-X(C))/RBC4
         VIXEN5= ZPP*ZE*VIXEN4 + (1.0D0+ZE*DMUZ4)*ZE*DMUX2 + XPP*(XE*VIXEN4-DMUZ4)
     1          + XE*DMUZ4*(XE*DMUX2-1.0D0-MU) + YPP*YE*VIXEN4 + YE**2*DMUZ4*DMUX2 - DQZ4*DQX2
         VIXEN3=(VIXEN1/VIXEN2)*VIXEN5 + DPX2*DQZ4
         VIXEN6= ZP*ZE*VIXEN4 + (1.0D0+ZE*DMUZ4)*ZE*DLAMBDAX2 + XP*(XE*VIXEN4-DMUZ4)
     1          + XE*DMUZ4*(XE*DLAMBDAX2-1.0D0-LAMBDA) + YP*YE*VIXEN4 + YE**2*DMUZ4*DLAMBDAX2
         DNX2= (V*VIXEN6 + DVX2*DUZ4 - U*VIXEN3 - DUX2*DVZ4)/V**2 - 2.0D0*DVX2*N/V

         IMPHELL(3*D,3*B-2)= IMPHELL(3*D,3*B-2)+ M*DNX2 + N*DMX2
         IMPHELL(3*B-2,3*D)= IMPHELL(3*B-2,3*D)+ M*DNX2 + N*DMX2
C
C DZ4DY2,DY2DZ4
C
         N=(V*DUZ4 - U*DVZ4)/V**2
         VIXEN4= -2.0D0*(Z(B)-Z(C))*(Y(B)-Y(C))/RBC4
         VIXEN5= ZPP*ZE*VIXEN4 + (1.0D0+ZE*DMUZ4)*ZE*DMUY2 + YPP*(YE*VIXEN4-DMUZ4)
     1          + YE*DMUZ4*(YE*DMUY2-1.0D0-MU) + XPP*XE*VIXEN4 + XE**2*DMUZ4*DMUY2 - DQZ4*DQY2
         VIXEN3=(VIXEN1/VIXEN2)*VIXEN5 + DPY2*DQZ4
         VIXEN6= ZP*ZE*VIXEN4 + (1.0D0+ZE*DMUZ4)*ZE*DLAMBDAY2 + YP*(YE*VIXEN4-DMUZ4)
     1          + YE*DMUZ4*(YE*DLAMBDAY2-1.0D0-LAMBDA) + XP*XE*VIXEN4 + XE**2*DMUZ4*DLAMBDAY2
         DNY2= (V*VIXEN6 + DVY2*DUZ4 - U*VIXEN3 - DUY2*DVZ4)/V**2 - 2.0D0*DVY2*N/V

         IMPHELL(3*D,3*B-1)= IMPHELL(3*D,3*B-1)+ M*DNY2 + N*DMY2
         IMPHELL(3*B-1,3*D)= IMPHELL(3*B-1,3*D)+ M*DNY2 + N*DMY2
C
C DX4DX3,DX3DX4
C
         N=(V*DUX4 - U*DVX4)/V**2
         VIXEN4= 2.0D0*(X(B)-X(C))**2/RBC4 - 1.0D0/RBC2
         VIXEN5= XPP*(XE*VIXEN4+DMUX4) + (1.0D0+XE*DMUX4)*(XE*DMUX3+MU) + YPP*YE*VIXEN4
     1          + YE**2*DMUX4*DMUX3 + ZPP*ZE*VIXEN4 + ZE**2*DMUX4*DMUX3 - DQX4*DQX3
         VIXEN3=(VIXEN1/VIXEN2)*VIXEN5 + DPX3*DQX4
         VIXEN6= XP*(XE*VIXEN4+DMUX4) + (1.0D0+XE*DMUX4)*(XE*DLAMBDAX3+LAMBDA) + VIXEN4*
     1          (YP*YE+ZP*ZE) + (YE**2+ZE**2)*DMUX4*DLAMBDAX3
         DNX3= (V*VIXEN6 + DVX3*DUX4 - U*VIXEN3 - DUX3*DVX4)/V**2 - 2.0D0*DVX3*N/V

         IMPHELL(3*D-2,3*C-2)= IMPHELL(3*D-2,3*C-2)+ M*DNX3 + N*DMX3
         IMPHELL(3*C-2,3*D-2)= IMPHELL(3*C-2,3*D-2)+ M*DNX3 + N*DMX3
C
C DY4DY3,DY3DY4
C
         N=(V*DUY4 - U*DVY4)/V**2
         VIXEN4= 2.0D0*(Y(B)-Y(C))**2/RBC4 - 1.0D0/RBC2
         VIXEN5= YPP*(YE*VIXEN4+DMUY4) + (1.0D0+YE*DMUY4)*(YE*DMUY3+MU) + XPP*XE*VIXEN4
     1          + XE**2*DMUY4*DMUY3 + ZPP*ZE*VIXEN4 + ZE**2*DMUY4*DMUY3 - DQY4*DQY3
         VIXEN3=(VIXEN1/VIXEN2)*VIXEN5 + DPY3*DQY4
         VIXEN6= YP*(YE*VIXEN4+DMUY4) + (1.0D0+YE*DMUY4)*(YE*DLAMBDAY3+LAMBDA) + VIXEN4*
     1          (XP*XE+ZP*ZE) + (XE**2+ZE**2)*DMUY4*DLAMBDAY3
         DNY3= (V*VIXEN6 + DVY3*DUY4 - U*VIXEN3 - DUY3*DVY4)/V**2 - 2.0D0*DVY3*N/V

         IMPHELL(3*D-1,3*C-1)= IMPHELL(3*D-1,3*C-1)+ M*DNY3 + N*DMY3
         IMPHELL(3*C-1,3*D-1)= IMPHELL(3*C-1,3*D-1)+ M*DNY3 + N*DMY3
C
C DZ4DZ3,DZ3DZ4
C
         N=(V*DUZ4 - U*DVZ4)/V**2
         VIXEN4= 2.0D0*(Z(B)-Z(C))**2/RBC4 - 1.0D0/RBC2
         VIXEN5= ZPP*(ZE*VIXEN4+DMUZ4) + (1.0D0+ZE*DMUZ4)*(ZE*DMUZ3+MU) + XPP*XE*VIXEN4
     1          + XE**2*DMUZ4*DMUZ3 + YPP*YE*VIXEN4 + YE**2*DMUZ4*DMUZ3 - DQZ4*DQZ3
         VIXEN3=(VIXEN1/VIXEN2)*VIXEN5 + DPZ3*DQZ4
         VIXEN6= ZP*(ZE*VIXEN4+DMUZ4) + (1.0D0+ZE*DMUZ4)*(ZE*DLAMBDAZ3+LAMBDA) + VIXEN4*
     1          (XP*XE+YP*YE) + (XE**2+YE**2)*DMUZ4*DLAMBDAZ3
         DNZ3= (V*VIXEN6 + DVZ3*DUZ4 - U*VIXEN3 - DUZ3*DVZ4)/V**2 - 2.0D0*DVZ3*N/V

         IMPHELL(3*D,3*C)= IMPHELL(3*D,3*C)+ M*DNZ3 + N*DMZ3
         IMPHELL(3*C,3*D)= IMPHELL(3*C,3*D)+ M*DNZ3 + N*DMZ3
C
C DX4DY3,DY3DX4
C
         N=(V*DUX4 - U*DVX4)/V**2
         VIXEN4= 2.0D0*(X(B)-X(C))*(Y(B)-Y(C))/RBC4
         VIXEN5= XPP*XE*VIXEN4 + (1.0D0+XE*DMUX4)*XE*DMUY3 + YPP*(YE*VIXEN4+DMUX4)
     1          + YE*DMUX4*(YE*DMUY3+MU) + ZPP*ZE*VIXEN4 + ZE**2*DMUX4*DMUY3 - DQX4*DQY3
         VIXEN3=(VIXEN1/VIXEN2)*VIXEN5 + DPY3*DQX4
         VIXEN6= XP*XE*VIXEN4 + (1.0D0+XE*DMUX4)*XE*DLAMBDAY3 + YP*(YE*VIXEN4+DMUX4)
     1          + YE*DMUX4*(YE*DLAMBDAY3+LAMBDA) + ZP*ZE*VIXEN4 + ZE**2*DMUX4*DLAMBDAY3
         DNY3= (V*VIXEN6 + DVY3*DUX4 - U*VIXEN3 - DUY3*DVX4)/V**2 - 2.0D0*DVY3*N/V

         IMPHELL(3*D-2,3*C-1)= IMPHELL(3*D-2,3*C-1) + M*DNY3 + N*DMY3
         IMPHELL(3*C-1,3*D-2)= IMPHELL(3*C-1,3*D-2) + M*DNY3 + N*DMY3
C
C DX4DZ3,DZ3DX4
C
         N=(V*DUX4 - U*DVX4)/V**2
         VIXEN4= 2.0D0*(X(B)-X(C))*(Z(B)-Z(C))/RBC4
         VIXEN5= XPP*XE*VIXEN4 + (1.0D0+XE*DMUX4)*XE*DMUZ3 + ZPP*(ZE*VIXEN4+DMUX4)
     1          + ZE*DMUX4*(ZE*DMUZ3+MU) + YPP*YE*VIXEN4 + YE**2*DMUX4*DMUZ3 - DQX4*DQZ3
         VIXEN3=(VIXEN1/VIXEN2)*VIXEN5 + DPZ3*DQX4
         VIXEN6= XP*XE*VIXEN4 + (1.0D0+XE*DMUX4)*XE*DLAMBDAZ3 + ZP*(ZE*VIXEN4+DMUX4)
     1          + ZE*DMUX4*(ZE*DLAMBDAZ3+LAMBDA) + YP*YE*VIXEN4 + YE**2*DMUX4*DLAMBDAZ3
         DNZ3= (V*VIXEN6 + DVZ3*DUX4 - U*VIXEN3 - DUZ3*DVX4)/V**2 - 2.0D0*DVZ3*N/V

         IMPHELL(3*D-2,3*C)= IMPHELL(3*D-2,3*C) + M*DNZ3 + N*DMZ3
         IMPHELL(3*C,3*D-2)= IMPHELL(3*C,3*D-2) + M*DNZ3 + N*DMZ3
C
C DY4DZ3,DZ3DY4
C
         N=(V*DUY4 - U*DVY4)/V**2
         VIXEN4= 2.0D0*(Y(B)-Y(C))*(Z(B)-Z(C))/RBC4
         VIXEN5= YPP*YE*VIXEN4 + (1.0D0+YE*DMUY4)*YE*DMUZ3 + ZPP*(ZE*VIXEN4+DMUY4)
     1          + ZE*DMUY4*(ZE*DMUZ3+MU) + XPP*XE*VIXEN4 + XE**2*DMUY4*DMUZ3 - DQY4*DQZ3
         VIXEN3=(VIXEN1/VIXEN2)*VIXEN5 + DPZ3*DQY4
         VIXEN6= YP*YE*VIXEN4 + (1.0D0+YE*DMUY4)*YE*DLAMBDAZ3 + ZP*(ZE*VIXEN4+DMUY4)
     1          + ZE*DMUY4*(ZE*DLAMBDAZ3+LAMBDA) + XP*XE*VIXEN4 + XE**2*DMUY4*DLAMBDAZ3
         DNZ3= (V*VIXEN6 + DVZ3*DUY4 - U*VIXEN3 - DUZ3*DVY4)/V**2 - 2.0D0*DVZ3*N/V

         IMPHELL(3*D-1,3*C)= IMPHELL(3*D-1,3*C) + M*DNZ3 + N*DMZ3
         IMPHELL(3*C,3*D-1)= IMPHELL(3*C,3*D-1) + M*DNZ3 + N*DMZ3
C
C DY4DX3,DX3DY4
C
         N=(V*DUY4 - U*DVY4)/V**2
         VIXEN4= 2.0D0*(Y(B)-Y(C))*(X(B)-X(C))/RBC4
         VIXEN5= YPP*YE*VIXEN4 + (1.0D0+YE*DMUY4)*YE*DMUX3 + XPP*(XE*VIXEN4+DMUY4)
     1          + XE*DMUY4*(XE*DMUX3+MU) + ZPP*ZE*VIXEN4 + ZE**2*DMUY4*DMUX3 - DQY4*DQX3
         VIXEN3=(VIXEN1/VIXEN2)*VIXEN5 + DPX3*DQY4
         VIXEN6= YP*YE*VIXEN4 + (1.0D0+YE*DMUY4)*YE*DLAMBDAX3 + XP*(XE*VIXEN4+DMUY4)
     1          + XE*DMUY4*(XE*DLAMBDAX3+LAMBDA) + ZP*ZE*VIXEN4 + ZE**2*DMUY4*DLAMBDAX3
         DNX3= (V*VIXEN6 + DVX3*DUY4 - U*VIXEN3 - DUX3*DVY4)/V**2 - 2.0D0*DVX3*N/V

         IMPHELL(3*D-1,3*C-2)= IMPHELL(3*D-1,3*C-2) + M*DNX3 + N*DMX3
         IMPHELL(3*C-2,3*D-1)= IMPHELL(3*C-2,3*D-1) + M*DNX3 + N*DMX3
C
C DZ4DX3,DX3DZ4
C
         N=(V*DUZ4 - U*DVZ4)/V**2
         VIXEN4= 2.0D0*(Z(B)-Z(C))*(X(B)-X(C))/RBC4
         VIXEN5= ZPP*ZE*VIXEN4 + (1.0D0+ZE*DMUZ4)*ZE*DMUX3 + XPP*(XE*VIXEN4+DMUZ4)
     1          + XE*DMUZ4*(XE*DMUX3+MU) + YPP*YE*VIXEN4 + YE**2*DMUZ4*DMUX3 - DQZ4*DQX3
         VIXEN3=(VIXEN1/VIXEN2)*VIXEN5 + DPX3*DQZ4
         VIXEN6= ZP*ZE*VIXEN4 + (1.0D0+ZE*DMUZ4)*ZE*DLAMBDAX3 + XP*(XE*VIXEN4+DMUZ4)
     1          + XE*DMUZ4*(XE*DLAMBDAX3+LAMBDA) + YP*YE*VIXEN4 + YE**2*DMUZ4*DLAMBDAX3
         DNX3= (V*VIXEN6 + DVX3*DUZ4 - U*VIXEN3 - DUX3*DVZ4)/V**2 - 2.0D0*DVX3*N/V

         IMPHELL(3*D,3*C-2)= IMPHELL(3*D,3*C-2) + M*DNX3 + N*DMX3
         IMPHELL(3*C-2,3*D)= IMPHELL(3*C-2,3*D) + M*DNX3 + N*DMX3
C
C DZ4DY3,DY3DZ4
C
         N=(V*DUZ4 - U*DVZ4)/V**2
         VIXEN4= 2.0D0*(Z(B)-Z(C))*(Y(B)-Y(C))/RBC4
         VIXEN5= ZPP*ZE*VIXEN4 + (1.0D0+ZE*DMUZ4)*ZE*DMUY3 + YPP*(YE*VIXEN4+DMUZ4)
     1          + YE*DMUZ4*(YE*DMUY3+MU) + XPP*XE*VIXEN4 + XE**2*DMUZ4*DMUY3 - DQZ4*DQY3
         VIXEN3=(VIXEN1/VIXEN2)*VIXEN5 + DPY3*DQZ4
         VIXEN6= ZP*ZE*VIXEN4 + (1.0D0+ZE*DMUZ4)*ZE*DLAMBDAY3 + YP*(YE*VIXEN4+DMUZ4)
     1          + YE*DMUZ4*(YE*DLAMBDAY3+LAMBDA) + XP*XE*VIXEN4 + XE**2*DMUZ4*DLAMBDAY3
         DNY3= (V*VIXEN6 + DVY3*DUZ4 - U*VIXEN3 - DUY3*DVZ4)/V**2 - 2.0D0*DVY3*N/V

         IMPHELL(3*D,3*C-1)= IMPHELL(3*D,3*C-1) + M*DNY3 + N*DMY3
         IMPHELL(3*C-1,3*D)= IMPHELL(3*C-1,3*D) + M*DNY3 + N*DMY3
C
C DX3DX3
C
         N=(V*DUX3 - U*DVX3)/V**2
         VIXEN4= 8.0D0*(XE**2-RBC2)*MUCRAP/RBC6 - 4.0D0*XE*(X(B)-X(D))/RBC4
         VIXEN5= 8.0D0*(XE**2-RBC2)*LAMBDACRAP/RBC6 - 4.0D0*XE*(X(B)-X(A))/RBC4
         VIXEN7= XPP*(2.0D0*DMUX3+XE*VIXEN4) + (MU+XE*DMUX3)**2 + VIXEN4*(YPP*YE+ZPP*ZE)
     1          + DMUX3**2*(YE**2+ZE**2) - DQX3**2
         VIXEN8= XP*(2.0D0*DLAMBDAX3+XE*VIXEN5) + (LAMBDA+XE*DLAMBDAX3)**2 + VIXEN5*(YP*YE+ZP*ZE)
     1          + DLAMBDAX3**2*(YE**2+ZE**2) - DPX3**2
         VIXEN3= (VIXEN1/VIXEN2)*VIXEN7 + (VIXEN2/VIXEN1)*VIXEN8 + 2.0D0*DPX3*DQX3
         VIXEN6= XP*(2.0D0*DMUX3+XE*VIXEN4) + 2.0D0*(MU+XE*DMUX3)*(LAMBDA+XE*DLAMBDAX3) 
     1          + XPP*(2.0D0*DLAMBDAX3+XE*VIXEN5) + VIXEN4*(YP*YE+ZP*ZE) 
     2          + VIXEN5*(YPP*YE+ZPP*ZE) + 2.0D0*DLAMBDAX3*DMUX3*(YE**2+ZE**2)
         DNX3= (V*VIXEN6 - U*VIXEN3)/V**2 - 2.0D0*DVX3*N/V

         IMPHELL(3*C-2,3*C-2)= IMPHELL(3*C-2,3*C-2)+ M*DNX3 + N*DMX3
C
C DY3DY3
C
         N=(V*DUY3 - U*DVY3)/V**2
         VIXEN4= 8.0D0*(XE**2-RBC2)*MUCRAP/RBC6 - 4.0D0*XE*(X(B)-X(D))/RBC4
         VIXEN5= 8.0D0*(XE**2-RBC2)*LAMBDACRAP/RBC6 - 4.0D0*XE*(X(B)-X(A))/RBC4
         VIXEN7= YPP*(2.0D0*DMUY3+YE*VIXEN4) + (MU+YE*DMUY3)**2 + VIXEN4*(XPP*XE+ZPP*ZE)
     1          + DMUY3**2*(XE**2+ZE**2) - DQY3**2
         VIXEN8= YP*(2.0D0*DLAMBDAY3+YE*VIXEN5) + (LAMBDA+YE*DLAMBDAY3)**2 + VIXEN5*(XP*XE+ZP*ZE)
     1          + DLAMBDAY3**2*(XE**2+ZE**2) - DPY3**2
         VIXEN3= (VIXEN1/VIXEN2)*VIXEN7 + (VIXEN2/VIXEN1)*VIXEN8 + 2.0D0*DPY3*DQY3
         VIXEN6= YP*(2.0D0*DMUY3+YE*VIXEN4) + 2.0D0*(MU+YE*DMUY3)*(LAMBDA+YE*DLAMBDAY3)
     1          + YPP*(2.0D0*DLAMBDAY3+YE*VIXEN5) + VIXEN4*(XP*XE+ZP*ZE)
     2          + VIXEN5*(XPP*XE+ZPP*ZE) + 2.0D0*DLAMBDAY3*DMUY3*(XE**2+ZE**2)
         DNY3= (V*VIXEN6 - U*VIXEN3)/V**2 - 2.0D0*DVY3*N/V

         IMPHELL(3*C-1,3*C-1)= IMPHELL(3*C-1,3*C-1)+ M*DNY3 + N*DMY3
C
C DZ3DZ3
C
         N=(V*DUZ3 - U*DVZ3)/V**2
         VIXEN4= 8.0D0*(XE**2-RBC2)*MUCRAP/RBC6 - 4.0D0*XE*(X(B)-X(D))/RBC4
         VIXEN5= 8.0D0*(XE**2-RBC2)*LAMBDACRAP/RBC6 - 4.0D0*XE*(X(B)-X(A))/RBC4
         VIXEN7= ZPP*(2.0D0*DMUZ3+ZE*VIXEN4) + (MU+ZE*DMUZ3)**2 + VIXEN4*(XPP*XE+YPP*YE)
     1          + DMUZ3**2*(XE**2+YE**2) - DQZ3**2
         VIXEN8= ZP*(2.0D0*DLAMBDAZ3+ZE*VIXEN5) + (LAMBDA+ZE*DLAMBDAZ3)**2 + VIXEN5*(XP*XE+YP*YE)
     1          + DLAMBDAZ3**2*(XE**2+YE**2) - DPZ3**2
         VIXEN3= (VIXEN1/VIXEN2)*VIXEN7 + (VIXEN2/VIXEN1)*VIXEN8 + 2.0D0*DPZ3*DQZ3
         VIXEN6= ZP*(2.0D0*DMUZ3+ZE*VIXEN4) + 2.0D0*(MU+ZE*DMUZ3)*(LAMBDA+ZE*DLAMBDAZ3)
     1          + ZPP*(2.0D0*DLAMBDAZ3+ZE*VIXEN5) + VIXEN4*(XP*XE+YP*YE)
     2          + VIXEN5*(XPP*XE+YPP*YE) + 2.0D0*DLAMBDAZ3*DMUZ3*(XE**2+YE**2)
         DNZ3= (V*VIXEN6 - U*VIXEN3)/V**2 - 2.0D0*DVZ3*N/V

         IMPHELL(3*C,3*C)= IMPHELL(3*C,3*C)+ M*DNZ3 + N*DMZ3
C
C DX3DY3,DY3DX3
C
         N=(V*DUX3 - U*DVX3)/V**2
         VIXEN4= 8.0D0*MUCRAP*XE*YE/RBC6 - 2.0D0*(YE*(X(B)-X(D))+XE*(Y(B)-Y(D)))/RBC4
         VIXEN5= 8.0D0*MUCRAP*XE*YE/RBC6 - 2.0D0*(YE*(X(B)-X(A))+XE*(Y(B)-Y(A)))/RBC4
         VIXEN7= XPP*(DMUY3+XE*VIXEN4) + (MU+XE*DMUX3)*XE*DMUY3 + YPP*(YE*VIXEN4+DMUX3)
     1          + YE*DMUX3*(MU+YE*DMUY3) + ZPP*ZE*VIXEN4 + ZE**2*DMUX3*DMUY3 - DQX3*DQY3
         VIXEN8= XP*(DLAMBDAY3+XE*VIXEN5) + (LAMBDA+XE*DLAMBDAX3)*XE*DLAMBDAY3 
     1          + YP*(YE*VIXEN5 + DLAMBDAX3) + YE*DLAMBDAX3*(LAMBDA+YE*DLAMBDAY3) 
     2          + ZP*ZE*VIXEN5 + ZE**2*DLAMBDAX3*DLAMBDAY3 - DPX3*DPY3
         VIXEN3= (VIXEN1/VIXEN2)*VIXEN7 + (VIXEN2/VIXEN1)*VIXEN8 + DPX3*DQY3 + DPY3*DQX3
         VIXEN6= XP*(DMUY3+XE*VIXEN4) + (MU+XE*DMUX3)*XE*DLAMBDAY3 + XPP*(DLAMBDAY3+XE*VIXEN5)
     1          + (LAMBDA+XE*DLAMBDAX3)*XE*DMUY3 + YP*(YE*VIXEN4+DMUX3) 
     2          + YE*DMUX3*(YE*DLAMBDAY3+LAMBDA) + YPP*(YE*VIXEN5 + DLAMBDAX3) 
     3          + YE*DLAMBDAX3*(MU+YE*DMUY3) + ZP*ZE*VIXEN4 + ZPP*ZE*VIXEN5
     4          + ZE**2*(DLAMBDAX3*DMUY3+DMUX3*DLAMBDAY3)
         DNY3= (V*VIXEN6 + DVY3*DUX3 - U*VIXEN3 - DUY3*DVX3)/V**2 - 2.0D0*DVY3*N/V

         IMPHELL(3*C-2,3*C-1)= IMPHELL(3*C-2,3*C-1)+ M*DNY3 + N*DMY3
         IMPHELL(3*C-1,3*C-2)= IMPHELL(3*C-1,3*C-2)+ M*DNY3 + N*DMY3
C
C DX3DZ3,DZ3DX3
C
         N=(V*DUX3 - U*DVX3)/V**2
         VIXEN4= 8.0D0*MUCRAP*XE*ZE/RBC6 - 2.0D0*(ZE*(X(B)-X(D))+XE*(Z(B)-Z(D)))/RBC4
         VIXEN5= 8.0D0*MUCRAP*XE*ZE/RBC6 - 2.0D0*(ZE*(X(B)-X(A))+XE*(Z(B)-Z(A)))/RBC4
         VIXEN7= XPP*(DMUZ3+XE*VIXEN4) + (MU+XE*DMUX3)*XE*DMUZ3 + ZPP*(ZE*VIXEN4+DMUX3)
     1          + ZE*DMUX3*(MU+ZE*DMUZ3) + YPP*YE*VIXEN4 + YE**2*DMUX3*DMUZ3 - DQX3*DQZ3
         VIXEN8= XP*(DLAMBDAZ3+XE*VIXEN5) + (LAMBDA+XE*DLAMBDAX3)*XE*DLAMBDAZ3
     1          + ZP*(ZE*VIXEN5 + DLAMBDAX3) + ZE*DLAMBDAX3*(LAMBDA+ZE*DLAMBDAZ3)
     2          + YP*YE*VIXEN5 + YE**2*DLAMBDAX3*DLAMBDAZ3 - DPX3*DPZ3
         VIXEN3= (VIXEN1/VIXEN2)*VIXEN7 + (VIXEN2/VIXEN1)*VIXEN8 + DPX3*DQZ3 + DPZ3*DQX3
         VIXEN6= XP*(DMUZ3+XE*VIXEN4) + (MU+XE*DMUX3)*XE*DLAMBDAZ3 + XPP*(DLAMBDAZ3+XE*VIXEN5)
     1          + (LAMBDA+XE*DLAMBDAX3)*XE*DMUZ3 + ZP*(ZE*VIXEN4+DMUX3)
     2          + ZE*DMUX3*(ZE*DLAMBDAZ3+LAMBDA) + ZPP*(ZE*VIXEN5 + DLAMBDAX3)
     3          + ZE*DLAMBDAX3*(MU+ZE*DMUZ3) + YP*YE*VIXEN4 + YPP*YE*VIXEN5
     4          + YE**2*(DLAMBDAX3*DMUZ3+DMUX3*DLAMBDAZ3)
         DNZ3= (V*VIXEN6 + DVZ3*DUX3 - U*VIXEN3 - DUZ3*DVX3)/V**2 - 2.0D0*DVZ3*N/V

         IMPHELL(3*C-2,3*C)= IMPHELL(3*C-2,3*C)+ M*DNZ3 + N*DMZ3
         IMPHELL(3*C,3*C-2)= IMPHELL(3*C,3*C-2)+ M*DNZ3 + N*DMZ3
C
C DY3DZ3,DZ3DY3
C
         N=(V*DUY3 - U*DVY3)/V**2
         VIXEN4= 8.0D0*MUCRAP*YE*ZE/RBC6 - 2.0D0*(ZE*(Y(B)-Y(D))+YE*(Z(B)-Z(D)))/RBC4
         VIXEN5= 8.0D0*MUCRAP*YE*ZE/RBC6 - 2.0D0*(ZE*(Y(B)-Y(A))+YE*(Z(B)-Z(A)))/RBC4
         VIXEN7= YPP*(DMUZ3+YE*VIXEN4) + (MU+YE*DMUY3)*YE*DMUZ3 + ZPP*(ZE*VIXEN4+DMUY3)
     1          + ZE*DMUY3*(MU+ZE*DMUZ3) + XPP*XE*VIXEN4 + XE**2*DMUY3*DMUZ3 - DQY3*DQZ3
         VIXEN8= YP*(DLAMBDAZ3+YE*VIXEN5) + (LAMBDA+YE*DLAMBDAY3)*YE*DLAMBDAZ3
     1          + ZP*(ZE*VIXEN5 + DLAMBDAY3) + ZE*DLAMBDAY3*(LAMBDA+ZE*DLAMBDAZ3)
     2          + XP*XE*VIXEN5 + XE**2*DLAMBDAY3*DLAMBDAZ3 - DPY3*DPZ3
         VIXEN3= (VIXEN1/VIXEN2)*VIXEN7 + (VIXEN2/VIXEN1)*VIXEN8 + DPY3*DQZ3 + DPZ3*DQY3
         VIXEN6= YP*(DMUZ3+YE*VIXEN4) + (MU+YE*DMUY3)*YE*DLAMBDAZ3 + YPP*(DLAMBDAZ3+YE*VIXEN5)
     1          + (LAMBDA+YE*DLAMBDAY3)*YE*DMUZ3 + ZP*(ZE*VIXEN4+DMUY3)
     2          + ZE*DMUY3*(ZE*DLAMBDAZ3+LAMBDA) + ZPP*(ZE*VIXEN5 + DLAMBDAY3)
     3          + ZE*DLAMBDAY3*(MU+ZE*DMUZ3) + XP*XE*VIXEN4 + XPP*XE*VIXEN5
     4          + XE**2*(DLAMBDAY3*DMUZ3+DMUY3*DLAMBDAZ3)
         DNZ3= (V*VIXEN6 + DVZ3*DUY3 - U*VIXEN3 - DUZ3*DVY3)/V**2 - 2.0D0*DVZ3*N/V

         IMPHELL(3*C-1,3*C)= IMPHELL(3*C-1,3*C)+ M*DNZ3 + N*DMZ3
         IMPHELL(3*C,3*C-1)= IMPHELL(3*C,3*C-1)+ M*DNZ3 + N*DMZ3
C
C DX2DX2
C
         N=(V*DUX2 - U*DVX2)/V**2
         VIXEN4= -2.0D0/RBC2 - (4.0D0*(X(B)-X(C))*(X(C)+X(D)+2.0D0*X(B))-
     1          2.0D0*MUCRAP)/RBC4 + 8.0D0*(X(B)-X(C))**2*MUCRAP/RBC6
         VIXEN5= -2.0D0/RBC2 - (4.0D0*(X(B)-X(C))*(X(C)+X(A)+2.0D0*X(B))-
     1    2.0D0*LAMBDACRAP)/RBC4 + 8.0D0*(X(B)-X(C))**2*LAMBDACRAP/RBC6
         VIXEN7= XP*(XE*VIXEN5-2.0D0*DLAMBDAX2) + (XE*DLAMBDAX2
     1    -1.0D0-LAMBDA)**2 + (YP*YE+ZP*ZE)*VIXEN5 - DPX2**2
     2    + (YE**2+ZE**2)*DLAMBDAX2**2
         VIXEN8= XPP*(XE*VIXEN4-2.0D0*DMUX2) + (XE*DMUX2-1.0D0-MU)
     1    **2 + (YPP*YE+ZPP*ZE)*VIXEN4 - DQX2**2
     2    + (YE**2+ZE**2)*DMUX2**2
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN7 + (VIXEN1/VIXEN2)*VIXEN8+2.0D0*DPX2*DQX2
         VIXEN6= XP*(XE*VIXEN4-2.0D0*DMUX2) + 2.0D0*(XE*DMUX2-1.0D0-MU)
     1    *(XE*DLAMBDAX2-1.0D0-LAMBDA) + XPP*(XE*VIXEN5-2.0D0*
     2    DLAMBDAX2) + VIXEN4*(YP*YE+ZP*ZE) + VIXEN5*
     3    (YPP*YE+ZPP*ZE) + 2.0D0*DMUX2*DLAMBDAX2*(YE
     4    **2+ZE**2)
         DNX2= (V*VIXEN6 - U*VIXEN3)/V**2 - 2.0D0*DVX2*N/V

         IMPHELL(3*B-2,3*B-2)= IMPHELL(3*B-2,3*B-2)+ M*DNX2 + N*DMX2
C
C DY2DY2
C
         N=(V*DUY2 - U*DVY2)/V**2
         VIXEN4= -2.0D0/RBC2 - (4.0D0*(Y(B)-Y(C))*(Y(C)+Y(D)+2.0D0*Y(B))-
     1          2.0D0*MUCRAP)/RBC4 + 8.0D0*(Y(B)-Y(C))**2*MUCRAP/RBC6
         VIXEN5= -2.0D0/RBC2 - (4.0D0*(Y(B)-Y(C))*(Y(C)+Y(A)+2.0D0*Y(B))-
     1    2.0D0*LAMBDACRAP)/RBC4 + 8.0D0*(Y(B)-Y(C))**2*LAMBDACRAP/RBC6
         VIXEN7= YP*(YE*VIXEN5-2.0D0*DLAMBDAY2) + (YE*DLAMBDAY2
     1    -1.0D0-LAMBDA)**2 + (XP*XE+ZP*ZE)*VIXEN5 - DPY2**2
     2    + (XE**2+ZE**2)*DLAMBDAY2**2
         VIXEN8= YPP*(YE*VIXEN4-2.0D0*DMUY2) + (YE*DMUY2-1.0D0-MU)
     1    **2 + (XPP*XE+ZPP*ZE)*VIXEN4 - DQY2**2
     2    + (XE**2+ZE**2)*DMUY2**2
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN7 + (VIXEN1/VIXEN2)*VIXEN8+2.0D0*DPY2*DQY2
         VIXEN6= YP*(YE*VIXEN4-2.0D0*DMUY2) + 2.0D0*(YE*DMUY2-1.0D0-MU)
     1    *(YE*DLAMBDAY2-1.0D0-LAMBDA) + YPP*(YE*VIXEN5-2.0D0*
     2    DLAMBDAY2) + VIXEN4*(XP*XE+ZP*ZE) + VIXEN5*
     3    (XPP*XE+ZPP*ZE) + 2.0D0*DMUY2*DLAMBDAY2*(XE
     4    **2+ZE**2)
         DNY2= (V*VIXEN6 - U*VIXEN3)/V**2 - 2.0D0*DVY2*N/V

         IMPHELL(3*B-1,3*B-1)= IMPHELL(3*B-1,3*B-1)+ M*DNY2 + N*DMY2
C
C DZ2DZ2
C
         N=(V*DUZ2 - U*DVZ2)/V**2
         VIXEN4= -2.0D0/RBC2 - (4.0D0*(Z(B)-Z(C))*(Z(C)+Z(D)+2.0D0*Z(B))-
     1          2.0D0*MUCRAP)/RBC4 + 8.0D0*(Z(B)-Z(C))**2*MUCRAP/RBC6
         VIXEN5= -2.0D0/RBC2 - (4.0D0*(Z(B)-Z(C))*(Z(C)+Z(A)+2.0D0*Z(B))-
     1    2.0D0*LAMBDACRAP)/RBC4 + 8.0D0*(Z(B)-Z(C))**2*LAMBDACRAP/RBC6
         VIXEN7= ZP*(ZE*VIXEN5-2.0D0*DLAMBDAZ2) + (ZE*DLAMBDAZ2
     1    -1.0D0-LAMBDA)**2 + (XP*XE+YP*YE)*VIXEN5 - DPZ2**2
     2    + (XE**2+YE**2)*DLAMBDAZ2**2
         VIXEN8= ZPP*(ZE*VIXEN4-2.0D0*DMUZ2) + (ZE*DMUZ2-1.0D0-MU)
     1    **2 + (XPP*XE+YPP*YE)*VIXEN4 - DQZ2**2
     2    + (XE**2+YE**2)*DMUZ2**2
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN7 + (VIXEN1/VIXEN2)*VIXEN8+2.0D0*DPZ2*DQZ2
         VIXEN6= ZP*(ZE*VIXEN4-2*DMUZ2) + 2.0D0*(ZE*DMUZ2-1.0D0-MU)
     1    *(ZE*DLAMBDAZ2-1.0D0-LAMBDA) + ZPP*(ZE*VIXEN5-2.0D0*
     2    DLAMBDAZ2) + VIXEN4*(XP*XE+YP*YE) + VIXEN5*
     3    (XPP*XE+YPP*YE) + 2.0D0*DMUZ2*DLAMBDAZ2*(XE
     4    **2+YE**2)
         DNZ2= (V*VIXEN6 - U*VIXEN3)/V**2 - 2.0D0*DVZ2*N/V

         IMPHELL(3*B,3*B)= IMPHELL(3*B,3*B)+ M*DNZ2 + N*DMZ2
C
C DX2DY2,DY2DX2
C
         N=(V*DUX2 - U*DVX2)/V**2
         VIXEN4= (-2.0D0*(Y(B)-Y(C))*(X(C)+X(D)-2.0D0*X(B))-2.0D0*(X(B)-X(C))*(Y(C)+Y(D)-2.0D0*Y(B)))/RBC4 
     1          + 8.0D0*(X(B)-X(C))*(Y(B)-Y(C))*MUCRAP
     2          /RBC6
         VIXEN5= (-2.0D0*(Y(B)-Y(C))*(X(C)+X(A)-2.0D0*X(B))-2.0D0*(X(B)-X(C))*(Y(C)+Y(A)-2.0D0*Y(B)))/RBC4 
     1          + 8.0D0*(X(B)-X(C))*(Y(B)-Y(C))
     2          *LAMBDACRAP/RBC6
         VIXEN7= XP*(XE*VIXEN5-DLAMBDAY2) + (XE*DLAMBDAX2-1.0D0-LAMBDA)*XE*DLAMBDAY2 + YP*(YE
     1          *VIXEN5-DLAMBDAX2) + YE*DLAMBDAX2*(YE*DLAMBDAY2-1.0D0-LAMBDA) + ZP*ZE*VIXEN5
     2          + ZE**2*DLAMBDAX2*DLAMBDAY2 - DPX2*DPY2
         VIXEN8= XPP*(XE*VIXEN4-DMUY2) + (XE*DMUX2-1.0D0-MU)*XE*DMUY2 + YPP*(YE
     1          *VIXEN4-DMUX2) + YE*DMUX2*(YE*DMUY2-1.0D0-MU) + ZPP*ZE*VIXEN4
     2          + ZE**2*DMUX2*DMUY2 - DQX2*DQY2
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN7 + (VIXEN1/VIXEN2)*VIXEN8 + DPX2*DQY2 + DQX2*DPY2
         VIXEN6= XP*(XE*VIXEN4-DMUY2) + (XE*DMUX2-1.0D0-MU)*XE*DLAMBDAY2 + XPP*(XE*VIXEN5-DLAMBDAY2)
     1          + (XE*DLAMBDAX2-1.0D0-LAMBDA)*XE*DMUY2 + YP*(YE*VIXEN4-DMUX2) + YE*DMUX2
     2          *(YE*DLAMBDAY2-1.0D0-LAMBDA) + YPP*(YE*VIXEN5-DLAMBDAX2) + YE*DLAMBDAX2
     3          *(YE*DMUY2-1.0D0-MU) + ZP*ZE*VIXEN4 + ZPP*ZE*VIXEN5 + ZE**2*(DMUX2*DLAMBDAY2
     4          +DMUY2*DLAMBDAX2)
         DNY2= (V*VIXEN6 + DVY2*DUX2 - U*VIXEN3 - DVX2*DUY2)/V**2 - 2.0D0*DVY2*N/V

         IMPHELL(3*B-2,3*B-1)= IMPHELL(3*B-2,3*B-1)+ M*DNY2 + N*DMY2
         IMPHELL(3*B-1,3*B-2)= IMPHELL(3*B-1,3*B-2)+ M*DNY2 + N*DMY2
C
C DX2DZ2,DZ2DX2
C
         N=(V*DUX2 - U*DVX2)/V**2
         VIXEN4= (-2.0D0*(Z(B)-Z(C))*(X(C)+X(D)-2.0D0*X(B))-2.0D0*(X(B)-X(C))*(Z(C)+Z(D)-2.0D0*Z(B)))/RBC4 
     1          + 8.0D0*(X(B)-X(C))*(Z(B)-Z(C))*MUCRAP
     2          /RBC6
         VIXEN5= (-2.0D0*(Z(B)-Z(C))*(X(C)+X(A)-2.0D0*X(B))-2.0D0*(X(B)-X(C))*(Z(C)+Z(A)-2.0D0*Z(B)))/RBC4 
     1          + 8.0D0*(X(B)-X(C))*(Z(B)-Z(C))
     2          *LAMBDACRAP/RBC6
         VIXEN7= XP*(XE*VIXEN5-DLAMBDAZ2) + (XE*DLAMBDAX2-1.0D0-LAMBDA)*XE*DLAMBDAZ2 + ZP*(ZE
     1          *VIXEN5-DLAMBDAX2) + ZE*DLAMBDAX2*(ZE*DLAMBDAZ2-1.0D0-LAMBDA) + YP*YE*VIXEN5
     2          + YE**2*DLAMBDAX2*DLAMBDAZ2 - DPX2*DPZ2
         VIXEN8= XPP*(XE*VIXEN4-DMUZ2) + (XE*DMUX2-1.0D0-MU)*XE*DMUZ2 + ZPP*(ZE
     1          *VIXEN4-DMUX2) + ZE*DMUX2*(ZE*DMUZ2-1.0D0-MU) + YPP*YE*VIXEN4
     2          + YE**2*DMUX2*DMUZ2 - DQX2*DQZ2
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN7 + (VIXEN1/VIXEN2)*VIXEN8 + DPX2*DQZ2 + DQX2*DPZ2
         VIXEN6= XP*(XE*VIXEN4-DMUZ2) + (XE*DMUX2-1.0D0-MU)*XE*DLAMBDAZ2 + XPP*(XE*VIXEN5-DLAMBDAZ2)
     1          + (XE*DLAMBDAX2-1.0D0-LAMBDA)*XE*DMUZ2 + ZP*(ZE*VIXEN4-DMUX2) + ZE*DMUX2
     2          *(ZE*DLAMBDAZ2-1.0D0-LAMBDA) + ZPP*(ZE*VIXEN5-DLAMBDAX2) + ZE*DLAMBDAX2
     3          *(ZE*DMUZ2-1.0D0-MU) + YP*YE*VIXEN4 + YPP*YE*VIXEN5 + YE**2*(DMUX2*DLAMBDAZ2
     4          +DMUZ2*DLAMBDAX2)
         DNZ2= (V*VIXEN6 + DVZ2*DUX2 - U*VIXEN3 - DVX2*DUZ2)/V**2 - 2.0D0*DVZ2*N/V

         IMPHELL(3*B-2,3*B)= IMPHELL(3*B-2,3*B)+ M*DNZ2 + N*DMZ2
         IMPHELL(3*B,3*B-2)= IMPHELL(3*B,3*B-2)+ M*DNZ2 + N*DMZ2
C
C DY2DZ2,DZ2DY2
C
         N=(V*DUY2 - U*DVY2)/V**2
         VIXEN4= (-2.0D0*(Z(B)-Z(C))*(Y(C)+Y(D)-2.0D0*Y(B))-2.0D0*(Y(B)-Y(C))*(Z(C)+Z(D)-2.0D0*Z(B)))/RBC4 
     1          + 8.0D0*(Y(B)-Y(C))*(Z(B)-Z(C))*MUCRAP
     2          /RBC6
         VIXEN5= (-2.0D0*(Z(B)-Z(C))*(Y(C)+Y(A)-2.0D0*Y(B))-2.0D0*(Y(B)-Y(C))*(Z(C)+Z(A)-2.0D0*Z(B)))/RBC4 
     1          + 8.0D0*(Y(B)-Y(C))*(Z(B)-Z(C))
     2          *LAMBDACRAP/RBC6
         VIXEN7= YP*(YE*VIXEN5-DLAMBDAZ2) + (YE*DLAMBDAY2-1.0D0-LAMBDA)*YE*DLAMBDAZ2 + ZP*(ZE
     1          *VIXEN5-DLAMBDAY2) + ZE*DLAMBDAY2*(ZE*DLAMBDAZ2-1.0D0-LAMBDA) + XP*XE*VIXEN5
     2          + XE**2*DLAMBDAY2*DLAMBDAZ2 - DPY2*DPZ2
         VIXEN8= YPP*(YE*VIXEN4-DMUZ2) + (YE*DMUY2-1.0D0-MU)*YE*DMUZ2 + ZPP*(ZE
     1          *VIXEN4-DMUY2) + ZE*DMUY2*(ZE*DMUZ2-1.0D0-MU) + XPP*XE*VIXEN4
     2          + XE**2*DMUY2*DMUZ2 - DQY2*DQZ2
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN7 + (VIXEN1/VIXEN2)*VIXEN8 + DPY2*DQZ2 + DQY2*DPZ2
         VIXEN6= YP*(YE*VIXEN4-DMUZ2) + (YE*DMUY2-1.0D0-MU)*YE*DLAMBDAZ2 + YPP*(YE*VIXEN5-DLAMBDAZ2)
     1          + (YE*DLAMBDAY2-1.0D0-LAMBDA)*YE*DMUZ2 + ZP*(ZE*VIXEN4-DMUY2) + ZE*DMUY2
     2          *(ZE*DLAMBDAZ2-1.0D0-LAMBDA) + ZPP*(ZE*VIXEN5-DLAMBDAY2) + ZE*DLAMBDAY2
     3          *(ZE*DMUZ2-1.0D0-MU) + XP*XE*VIXEN4 + XPP*XE*VIXEN5 + XE**2*(DMUY2*DLAMBDAZ2
     4          +DMUZ2*DLAMBDAY2)
         DNZ2= (V*VIXEN6 + DVZ2*DUY2 - U*VIXEN3 - DVY2*DUZ2)/V**2 - 2.0D0*DVZ2*N/V

         IMPHELL(3*B-1,3*B)= IMPHELL(3*B-1,3*B)+ M*DNZ2 + N*DMZ2
         IMPHELL(3*B,3*B-1)= IMPHELL(3*B,3*B-1)+ M*DNZ2 + N*DMZ2
C
C DX2DX3,DX3DX2
C
         N=(V*DUX2 - U*DVX2)/V**2
         VIXEN4= 1.0D0/RBC2 - (2.0D0*XE*(X(C)+X(D)-2.0D0*X(B))+2.0D0*((X(B)-X(C))*(X(B)-X(D))-MUCRAP))/RBC4 
     1          - 8.0D0*(X(B)-X(C))**2*MUCRAP/RBC6
         VIXEN5= 1.0D0/RBC2 - (2.0D0*XE*(X(C)+X(A)-2.0D0*X(B))+2.0D0*((X(B)-X(C))*(X(B)-X(A))-LAMBDACRAP))/RBC4 
     1          - 8.0D0*(X(B)-X(C))**2*LAMBDACRAP/RBC6
         VIXEN7= XP*(XE*VIXEN5+DLAMBDAX2-DLAMBDAX3) + (XE*DLAMBDAX2-1.0D0-LAMBDA)*(XE*DLAMBDAX3+LAMBDA)
     1          + (YP*YE+ZP*ZE)*VIXEN5 + (YE**2+ZE**2)*DLAMBDAX2*DLAMBDAX3 - DPX2*DPX3
         VIXEN8= XPP*(XE*VIXEN4+DMUX2-DMUX3) + (XE*DMUX2-1.0D0-MU)*(XE*DMUX3+MU)
     1          + (YPP*YE+ZPP*ZE)*VIXEN4 + (YE**2+ZE**2)*DMUX2*DMUX3 - DQX2*DQX3
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN7 + (VIXEN1/VIXEN2)*VIXEN8 + DPX2*DQX3 + DQX2*DPX3
         VIXEN6= XP*(XE*VIXEN4+DMUX2-DMUX3) + (XE*DMUX2-1.0D0-MU)*(XE*DLAMBDAX3+LAMBDA)
     1          + XPP*(XE*VIXEN5+DLAMBDAX2-DLAMBDAX3) + (XE*DLAMBDAX2-1.0D0-LAMBDA)*(XE*DMUX3+MU)
     2          + VIXEN4*(YP*YE+ZP*ZE) + VIXEN5*(YPP*YE+ZPP*ZE) 
     3          + (DMUX2*DLAMBDAX3+DMUX3*DLAMBDAX2)*(YE**2+ZE**2)
         DNX3= (V*VIXEN6 + DVX3*DUX2 - U*VIXEN3 - DUX3*DVX2)/V**2 - 2.0D0*DVX3*N/V

         IMPHELL(3*B-2,3*C-2)= IMPHELL(3*B-2,3*C-2)+ M*DNX3 + N*DMX3
         IMPHELL(3*C-2,3*B-2)= IMPHELL(3*C-2,3*B-2)+ M*DNX3 + N*DMX3
C
C DY2DY3,DY3DY2
C
         N=(V*DUY2 - U*DVY2)/V**2
         VIXEN4= 1.0D0/RBC2 - (2.0D0*YE*(Y(C)+Y(D)-2.0D0*Y(B))+2.0D0*((Y(B)-Y(C))*(Y(B)-Y(D))-MUCRAP))/RBC4 
     1          - 8.0D0*(Y(B)-Y(C))**2*MUCRAP/RBC6
         VIXEN5= 1.0D0/RBC2 - (2.0D0*YE*(Y(C)+Y(A)-2.0D0*Y(B))+2.0D0*((Y(B)-Y(C))*(Y(B)-Y(A))-LAMBDACRAP))/RBC4 
     1          - 8.0D0*(Y(B)-Y(C))**2*LAMBDACRAP/RBC6
         VIXEN7= YP*(YE*VIXEN5+DLAMBDAY2-DLAMBDAY3) + (YE*DLAMBDAY2-1.0D0-LAMBDA)*(YE*DLAMBDAY3+LAMBDA)
     1          + (XP*XE+ZP*ZE)*VIXEN5 + (XE**2+ZE**2)*DLAMBDAY2*DLAMBDAY3 - DPY2*DPY3
         VIXEN8= YPP*(YE*VIXEN4+DMUY2-DMUY3) + (YE*DMUY2-1.0D0-MU)*(YE*DMUY3+MU)
     1          + (XPP*XE+ZPP*ZE)*VIXEN4 + (XE**2+ZE**2)*DMUY2*DMUY3 - DQY2*DQY3
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN7 + (VIXEN1/VIXEN2)*VIXEN8 + DPY2*DQY3 + DQY2*DPY3
         VIXEN6= YP*(YE*VIXEN4+DMUY2-DMUY3) + (YE*DMUY2-1.0D0-MU)*(YE*DLAMBDAY3+LAMBDA)
     1          + YPP*(YE*VIXEN5+DLAMBDAY2-DLAMBDAY3) + (YE*DLAMBDAY2-1.0D0-LAMBDA)*(YE*DMUY3+MU)
     2          + VIXEN4*(XP*XE+ZP*ZE) + VIXEN5*(XPP*XE+ZPP*ZE) 
     3          + (DMUY2*DLAMBDAY3+DMUY3*DLAMBDAY2)*(XE**2+ZE**2)
         DNY3= (V*VIXEN6 + DVY3*DUY2 - U*VIXEN3 - DUY3*DVY2)/V**2 - 2.0D0*DVY3*N/V

         IMPHELL(3*B-1,3*C-1)= IMPHELL(3*B-1,3*C-1)+ M*DNY3 + N*DMY3
         IMPHELL(3*C-1,3*B-1)= IMPHELL(3*C-1,3*B-1)+ M*DNY3 + N*DMY3
C
C DZ2DZ3,DZ3DZ2
C
         N=(V*DUZ2 - U*DVZ2)/V**2
         VIXEN4= 1.0D0/RBC2 - (2.0D0*ZE*(Z(C)+Z(D)-2.0D0*Z(B))+2.0D0*((Z(B)-Z(C))*(Z(B)-Z(D))-MUCRAP))/RBC4 
     1          - 8.0D0*(Z(B)-Z(C))**2*MUCRAP/RBC6
         VIXEN5= 1.0D0/RBC2 - (2.0D0*ZE*(Z(C)+Z(A)-2.0D0*Z(B))+2.0D0*((Z(B)-Z(C))*(Z(B)-Z(A))-LAMBDACRAP))/RBC4 
     1          - 8.0D0*(Z(B)-Z(C))**2*LAMBDACRAP/RBC6
         VIXEN7= ZP*(ZE*VIXEN5+DLAMBDAZ2-DLAMBDAZ3) + (ZE*DLAMBDAZ2-1.0D0-LAMBDA)*(ZE*DLAMBDAZ3+LAMBDA)
     1          + (XP*XE+YP*YE)*VIXEN5 + (XE**2+YE**2)*DLAMBDAZ2*DLAMBDAZ3 - DPZ2*DPZ3
         VIXEN8= ZPP*(ZE*VIXEN4+DMUZ2-DMUZ3) + (ZE*DMUZ2-1.0D0-MU)*(ZE*DMUZ3+MU)
     1          + (XPP*XE+YPP*YE)*VIXEN4 + (XE**2+YE**2)*DMUZ2*DMUZ3 - DQZ2*DQZ3
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN7 + (VIXEN1/VIXEN2)*VIXEN8 + DPZ2*DQZ3 + DQZ2*DPZ3
         VIXEN6= ZP*(ZE*VIXEN4+DMUZ2-DMUZ3) + (ZE*DMUZ2-1.0D0-MU)*(ZE*DLAMBDAZ3+LAMBDA)
     1          + ZPP*(ZE*VIXEN5+DLAMBDAZ2-DLAMBDAZ3) + (ZE*DLAMBDAZ2-1.0D0-LAMBDA)*(ZE*DMUZ3+MU)
     2          + VIXEN4*(XP*XE+YP*YE) + VIXEN5*(XPP*XE+YPP*YE) 
     3          + (DMUZ2*DLAMBDAZ3+DMUZ3*DLAMBDAZ2)*(XE**2+YE**2)
         DNZ3= (V*VIXEN6 + DVZ3*DUZ2 - U*VIXEN3 - DUZ3*DVZ2)/V**2 - 2.0D0*DVZ3*N/V

         IMPHELL(3*B,3*C)= IMPHELL(3*B,3*C)+ M*DNZ3 + N*DMZ3
         IMPHELL(3*C,3*B)= IMPHELL(3*C,3*B)+ M*DNZ3 + N*DMZ3
C
C DX2DY3,DY3DX2
C
         N=(V*DUX2 - U*DVX2)/V**2
         VIXEN4= 8.0D0*(X(B)-X(C))*YE*MUCRAP/RBC6 - 2.0D0*((X(C)+X(D)-2.0D0*X(B))*YE
     1          +(X(B)-X(C))*(Y(B)-Y(D)))/RBC4
         VIXEN5= 8.0D0*(X(B)-X(C))*YE*LAMBDACRAP/RBC6 - 2.0D0*((X(C)+X(A)-2.0D0*X(B))*YE
     1          +(X(B)-X(C))*(Y(B)-Y(A)))/RBC4
         VIXEN7= XP*(XE*VIXEN5-DLAMBDAY3) + (XE*DLAMBDAX2-1.0D0-LAMBDA)*XE*DLAMBDAY3 + 
     1          YP*(YE*VIXEN5+DLAMBDAX2) + YE*DLAMBDAX2*(YE*DLAMBDAY3+LAMBDA) 
     2          + ZP*ZE*VIXEN5 + ZE**2*DLAMBDAX2*DLAMBDAY3 - DPX2*DPY3
         VIXEN8= XPP*(XE*VIXEN4-DMUY3) + (XE*DMUX2-1.0D0-MU)*XE*DMUY3 + 
     1          YPP*(YE*VIXEN4+DMUX2) + YE*DMUX2*(YE*DMUY3+MU) + 
     2          ZPP*ZE*VIXEN4 + ZE**2*DMUX2*DMUY3 - DQX2*DQY3
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN7 + (VIXEN1/VIXEN2)*VIXEN8 + DPX2*DQY3 + DQX2*DPY3
         VIXEN6= XP*(XE*VIXEN4-DMUY3) + (XE*DMUX2-1.0D0-MU)*XE*DLAMBDAY3
     1          + XPP*(XE*VIXEN5-DLAMBDAY3) + (XE*DLAMBDAX2-1.0D0-LAMBDA)*XE*DMUY3
     2          + YP*(YE*VIXEN4+DMUX2) + YE*DMUX2*(YE*DLAMBDAY3+LAMBDA)
     3          + YPP*(YE*VIXEN5+DLAMBDAX2) + YE*DLAMBDAX2*(YE*DMUY3+MU)
     4          + ZP*ZE*VIXEN4 + ZE**2*(DMUX2*DLAMBDAY3+DLAMBDAX2*DMUY3) + ZPP*ZE*VIXEN5
         DNY3= (V*VIXEN6 + DVY3*DUX2 - U*VIXEN3 - DUY3*DVX2)/V**2 - 2.0D0*DVY3*N/V

         IMPHELL(3*B-2,3*C-1)= IMPHELL(3*B-2,3*C-1)+ M*DNY3 + N*DMY3
         IMPHELL(3*C-1,3*B-2)= IMPHELL(3*C-1,3*B-2)+ M*DNY3 + N*DMY3
C
C DX2DZ3,DZ3DX2
C
         N=(V*DUX2 - U*DVX2)/V**2
         VIXEN4= 8.0D0*(X(B)-X(C))*ZE*MUCRAP/RBC6 - 2.0D0*((X(C)+X(D)-2.0D0*X(B))*ZE
     1          +(X(B)-X(C))*(Z(B)-Z(D)))/RBC4
         VIXEN5= 8.0D0*(X(B)-X(C))*ZE*LAMBDACRAP/RBC6 - 2.0D0*((X(C)+X(A)-2.0D0*X(B))*ZE
     1          +(X(B)-X(C))*(Z(B)-Z(A)))/RBC4
         VIXEN7= XP*(XE*VIXEN5-DLAMBDAZ3) + (XE*DLAMBDAX2-1.0D0-LAMBDA)*XE*DLAMBDAZ3 + 
     1          ZP*(ZE*VIXEN5+DLAMBDAX2) + ZE*DLAMBDAX2*(ZE*DLAMBDAZ3+LAMBDA) 
     2          + YP*YE*VIXEN5 + YE**2*DLAMBDAX2*DLAMBDAZ3 - DPX2*DPZ3
         VIXEN8= XPP*(XE*VIXEN4-DMUZ3) + (XE*DMUX2-1.0D0-MU)*XE*DMUZ3 + 
     1          ZPP*(ZE*VIXEN4+DMUX2) + ZE*DMUX2*(ZE*DMUZ3+MU) + 
     2          YPP*YE*VIXEN4 + YE**2*DMUX2*DMUZ3 - DQX2*DQZ3
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN7 + (VIXEN1/VIXEN2)*VIXEN8 + DPX2*DQZ3 + DQX2*DPZ3
         VIXEN6= XP*(XE*VIXEN4-DMUZ3) + (XE*DMUX2-1.0D0-MU)*XE*DLAMBDAZ3
     1          + XPP*(XE*VIXEN5-DLAMBDAZ3) + (XE*DLAMBDAX2-1.0D0-LAMBDA)*XE*DMUZ3
     2          + ZP*(ZE*VIXEN4+DMUX2) + ZE*DMUX2*(ZE*DLAMBDAZ3+LAMBDA)
     3          + ZPP*(ZE*VIXEN5+DLAMBDAX2) + ZE*DLAMBDAX2*(ZE*DMUZ3+MU)
     4          + YP*YE*VIXEN4 + YE**2*(DMUX2*DLAMBDAZ3+DLAMBDAX2*DMUZ3) + YPP*YE*VIXEN5
         DNZ3= (V*VIXEN6 + DVZ3*DUX2 - U*VIXEN3 - DUZ3*DVX2)/V**2 - 2.0D0*DVZ3*N/V

         IMPHELL(3*B-2,3*C)= IMPHELL(3*B-2,3*C)+ M*DNZ3 + N*DMZ3
         IMPHELL(3*C,3*B-2)= IMPHELL(3*C,3*B-2)+ M*DNZ3 + N*DMZ3
C
C DY2DZ3,DZ3DY2
C
         N=(V*DUY2 - U*DVY2)/V**2
         VIXEN4= 8.0D0*(Y(B)-Y(C))*ZE*MUCRAP/RBC6 - 2.0D0*((Y(C)+Y(D)-2.0D0*Y(B))*ZE
     1          +(Y(B)-Y(C))*(Z(B)-Z(D)))/RBC4
         VIXEN5= 8.0D0*(Y(B)-Y(C))*ZE*LAMBDACRAP/RBC6 - 2.0D0*((Y(C)+Y(A)-2.0D0*Y(B))*ZE
     1          +(Y(B)-Y(C))*(Z(B)-Z(A)))/RBC4
         VIXEN7= YP*(YE*VIXEN5-DLAMBDAZ3) + (YE*DLAMBDAY2-1.0D0-LAMBDA)*YE*DLAMBDAZ3 +
     1          ZP*(ZE*VIXEN5+DLAMBDAY2) + ZE*DLAMBDAY2*(ZE*DLAMBDAZ3+LAMBDA)
     2          + XP*XE*VIXEN5 + XE**2*DLAMBDAY2*DLAMBDAZ3 - DPY2*DPZ3
         VIXEN8= YPP*(YE*VIXEN4-DMUZ3) + (YE*DMUY2-1.0D0-MU)*YE*DMUZ3 +
     1          ZPP*(ZE*VIXEN4+DMUY2) + ZE*DMUY2*(ZE*DMUZ3+MU) +
     2          XPP*XE*VIXEN4 + XE**2*DMUY2*DMUZ3 - DQY2*DQZ3
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN7 + (VIXEN1/VIXEN2)*VIXEN8 + DPY2*DQZ3 + DQY2*DPZ3
         VIXEN6= YP*(YE*VIXEN4-DMUZ3) + (YE*DMUY2-1.0D0-MU)*YE*DLAMBDAZ3
     1          + YPP*(YE*VIXEN5-DLAMBDAZ3) + (YE*DLAMBDAY2-1.0D0-LAMBDA)*YE*DMUZ3
     2          + ZP*(ZE*VIXEN4+DMUY2) + ZE*DMUY2*(ZE*DLAMBDAZ3+LAMBDA)
     3          + ZPP*(ZE*VIXEN5+DLAMBDAY2) + ZE*DLAMBDAY2*(ZE*DMUZ3+MU)
     4          + XP*XE*VIXEN4 + XE**2*(DMUY2*DLAMBDAZ3+DLAMBDAY2*DMUZ3) + XPP*XE*VIXEN5
         DNZ3= (V*VIXEN6 + DVZ3*DUY2 - U*VIXEN3 - DUZ3*DVY2)/V**2 - 2.0D0*DVZ3*N/V

         IMPHELL(3*B-1,3*C)= IMPHELL(3*B-1,3*C)+ M*DNZ3 + N*DMZ3
         IMPHELL(3*C,3*B-1)= IMPHELL(3*C,3*B-1)+ M*DNZ3 + N*DMZ3
C
C DY2DX3,DX3DY2
C
         N=(V*DUY2 - U*DVY2)/V**2
         VIXEN4= 8.0D0*(Y(B)-Y(C))*XE*MUCRAP/RBC6 - 2.0D0*((Y(C)+Y(D)-2.0D0*Y(B))*XE
     1          +(Y(B)-Y(C))*(X(B)-X(D)))/RBC4
         VIXEN5= 8.0D0*(Y(B)-Y(C))*XE*LAMBDACRAP/RBC6 - 2.0D0*((Y(C)+Y(A)-2.0D0*Y(B))*XE
     1          +(Y(B)-Y(C))*(X(B)-X(A)))/RBC4
         VIXEN7= YP*(YE*VIXEN5-DLAMBDAX3) + (YE*DLAMBDAY2-1.0D0-LAMBDA)*YE*DLAMBDAX3 +
     1          XP*(XE*VIXEN5+DLAMBDAY2) + XE*DLAMBDAY2*(XE*DLAMBDAX3+LAMBDA)
     2          + ZP*ZE*VIXEN5 + ZE**2*DLAMBDAY2*DLAMBDAX3 - DPY2*DPX3
         VIXEN8= YPP*(YE*VIXEN4-DMUX3) + (YE*DMUY2-1.0D0-MU)*YE*DMUX3 +
     1          XPP*(XE*VIXEN4+DMUY2) + XE*DMUY2*(XE*DMUX3+MU) +
     2          ZPP*ZE*VIXEN4 + ZE**2*DMUY2*DMUX3 - DQY2*DQX3
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN7 + (VIXEN1/VIXEN2)*VIXEN8 + DPY2*DQX3 + DQY2*DPX3
         VIXEN6= YP*(YE*VIXEN4-DMUX3) + (YE*DMUY2-1.0D0-MU)*YE*DLAMBDAX3
     1          + YPP*(YE*VIXEN5-DLAMBDAX3) + (YE*DLAMBDAY2-1.0D0-LAMBDA)*YE*DMUX3
     2          + XP*(XE*VIXEN4+DMUY2) + XE*DMUY2*(XE*DLAMBDAX3+LAMBDA)
     3          + XPP*(XE*VIXEN5+DLAMBDAY2) + XE*DLAMBDAY2*(XE*DMUX3+MU)
     4          + ZP*ZE*VIXEN4 + ZE**2*(DMUY2*DLAMBDAX3+DLAMBDAY2*DMUX3) + ZPP*ZE*VIXEN5
         DNX3= (V*VIXEN6 + DVX3*DUY2 - U*VIXEN3 - DUX3*DVY2)/V**2 - 2.0D0*DVX3*N/V

         IMPHELL(3*B-1,3*C-2)= IMPHELL(3*B-1,3*C-2)+ M*DNX3 + N*DMX3
         IMPHELL(3*C-2,3*B-1)= IMPHELL(3*C-2,3*B-1)+ M*DNX3 + N*DMX3
C
C DZ2DX3,DX3DZ2
C
         N=(V*DUZ2 - U*DVZ2)/V**2
         VIXEN4= 8.0D0*(Z(B)-Z(C))*XE*MUCRAP/RBC6 - 2.0D0*((Z(C)+Z(D)-2.0D0*Z(B))*XE
     1          +(Z(B)-Z(C))*(X(B)-X(D)))/RBC4
         VIXEN5= 8.0D0*(Z(B)-Z(C))*XE*LAMBDACRAP/RBC6 - 2.0D0*((Z(C)+Z(A)-2.0D0*Z(B))*XE
     1          +(Z(B)-Z(C))*(X(B)-X(A)))/RBC4
         VIXEN7= ZP*(ZE*VIXEN5-DLAMBDAX3) + (ZE*DLAMBDAZ2-1.0D0-LAMBDA)*ZE*DLAMBDAX3 +
     1          XP*(XE*VIXEN5+DLAMBDAZ2) + XE*DLAMBDAZ2*(XE*DLAMBDAX3+LAMBDA)
     2          + YP*YE*VIXEN5 + YE**2*DLAMBDAZ2*DLAMBDAX3 - DPZ2*DPX3
         VIXEN8= ZPP*(ZE*VIXEN4-DMUX3) + (ZE*DMUZ2-1.0D0-MU)*ZE*DMUX3 +
     1          XPP*(XE*VIXEN4+DMUZ2) + XE*DMUZ2*(XE*DMUX3+MU) +
     2          YPP*YE*VIXEN4 + YE**2*DMUZ2*DMUX3 - DQZ2*DQX3
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN7 + (VIXEN1/VIXEN2)*VIXEN8 + DPZ2*DQX3 + DQZ2*DPX3
         VIXEN6= ZP*(ZE*VIXEN4-DMUX3) + (ZE*DMUZ2-1.0D0-MU)*ZE*DLAMBDAX3
     1          + ZPP*(ZE*VIXEN5-DLAMBDAX3) + (ZE*DLAMBDAZ2-1.0D0-LAMBDA)*ZE*DMUX3
     2          + XP*(XE*VIXEN4+DMUZ2) + XE*DMUZ2*(XE*DLAMBDAX3+LAMBDA)
     3          + XPP*(XE*VIXEN5+DLAMBDAZ2) + XE*DLAMBDAZ2*(XE*DMUX3+MU)
     4          + YP*YE*VIXEN4 + YE**2*(DMUZ2*DLAMBDAX3+DLAMBDAZ2*DMUX3) + YPP*YE*VIXEN5
         DNX3= (V*VIXEN6 + DVX3*DUZ2 - U*VIXEN3 - DUX3*DVZ2)/V**2 - 2.0D0*DVX3*N/V

         IMPHELL(3*B,3*C-2)= IMPHELL(3*B,3*C-2)+ M*DNX3 + N*DMX3
         IMPHELL(3*C-2,3*B)= IMPHELL(3*C-2,3*B)+ M*DNX3 + N*DMX3
C
C DZ2DY3,DY3DZ2
C
         N=(V*DUZ2 - U*DVZ2)/V**2
         VIXEN4= 8.0D0*(Z(B)-Z(C))*YE*MUCRAP/RBC6 - 2.0D0*((Z(C)+Z(D)-2.0D0*Z(B))*YE
     1          +(Z(B)-Z(C))*(Y(B)-Y(D)))/RBC4
         VIXEN5= 8.0D0*(Z(B)-Z(C))*YE*LAMBDACRAP/RBC6 - 2.0D0*((Z(C)+Z(A)-2.0D0*Z(B))*YE
     1          +(Z(B)-Z(C))*(Y(B)-Y(A)))/RBC4
         VIXEN7= ZP*(ZE*VIXEN5-DLAMBDAY3) + (ZE*DLAMBDAZ2-1.0D0-LAMBDA)*ZE*DLAMBDAY3 +
     1          YP*(YE*VIXEN5+DLAMBDAZ2) + YE*DLAMBDAZ2*(YE*DLAMBDAY3+LAMBDA)
     2          + XP*XE*VIXEN5 + XE**2*DLAMBDAZ2*DLAMBDAY3 - DPZ2*DPY3
         VIXEN8= ZPP*(ZE*VIXEN4-DMUY3) + (ZE*DMUZ2-1.0D0-MU)*ZE*DMUY3 +
     1          YPP*(YE*VIXEN4+DMUZ2) + YE*DMUZ2*(YE*DMUY3+MU) +
     2          XPP*XE*VIXEN4 + XE**2*DMUZ2*DMUY3 - DQZ2*DQY3
         VIXEN3= (VIXEN2/VIXEN1)*VIXEN7 + (VIXEN1/VIXEN2)*VIXEN8 + DPZ2*DQY3 + DQZ2*DPY3
         VIXEN6= ZP*(ZE*VIXEN4-DMUY3) + (ZE*DMUZ2-1.0D0-MU)*ZE*DLAMBDAY3
     1          + ZPP*(ZE*VIXEN5-DLAMBDAY3) + (ZE*DLAMBDAZ2-1.0D0-LAMBDA)*ZE*DMUY3
     2          + YP*(YE*VIXEN4+DMUZ2) + YE*DMUZ2*(YE*DLAMBDAY3+LAMBDA)
     3          + YPP*(YE*VIXEN5+DLAMBDAZ2) + YE*DLAMBDAZ2*(YE*DMUY3+MU)
     4          + XP*XE*VIXEN4 + XE**2*(DMUZ2*DLAMBDAY3+DLAMBDAZ2*DMUY3) + XPP*XE*VIXEN5
         DNY3= (V*VIXEN6 + DVY3*DUZ2 - U*VIXEN3 - DUY3*DVZ2)/V**2 - 2.0D0*DVY3*N/V

         IMPHELL(3*B,3*C-1)= IMPHELL(3*B,3*C-1)+ M*DNY3 + N*DMY3
         IMPHELL(3*C-1,3*B)= IMPHELL(3*C-1,3*B)+ M*DNY3 + N*DMY3
      END DO

      DO I=1,3*ATOMS
         DO J=1,3*ATOMS
            HELL(I,J)=BONDHELL(I,J)+ANGLEHELL(I,J)+TORSHELL(I,J)+IMPHELL(I,J)+QHELL(I,J)+VDWHELL(I,J)
         END DO
      END DO

      RETURN

      END


      SUBROUTINE HAIRYHELL(U,V)
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE
      DOUBLE PRECISION      DUX1,DUX2,DUX3,DUX4,DUY1,DUY2,DUY3,DUY4,DUZ1,DUZ2,DUZ3,DUZ4
      COMMON /DU/           DUX1,DUX2,DUX3,DUX4,DUY1,DUY2,DUY3,DUY4,DUZ1,DUZ2,DUZ3,DUZ4
      DOUBLE PRECISION      DVX1,DVX2,DVX3,DVX4,DVY1,DVY2,DVY3,DVY4,DVZ1,DVZ2,DVZ3,DVZ4
      COMMON /DV/           DVX1,DVX2,DVX3,DVX4,DVY1,DVY2,DVY3,DVY4,DVZ1,DVZ2,DVZ3,DVZ4
      DOUBLE PRECISION      DMX1,DMX2,DMX3,DMX4,DMY1,DMY2,DMY3,DMY4,DMZ1,DMZ2,DMZ3,DMZ4
      COMMON /DM/           DMX1,DMX2,DMX3,DMX4,DMY1,DMY2,DMY3,DMY4,DMZ1,DMZ2,DMZ3,DMZ4
      DOUBLE PRECISION      U,V
      DOUBLE PRECISION      M,N
      COMMON /MN/           M,N

C      PRINT *,'IN HAIRYHELL',I
      M=0.0D0
      DMX1=0.0D0
      DMX2=0.0D0
      DMX3=0.0D0
      DMX4=0.0D0
      DMY1=0.0D0
      DMY2=0.0D0
      DMY3=0.0D0
      DMY4=0.0D0
      DMZ1=0.0D0
      DMZ2=0.0D0
      DMZ3=0.0D0
      DMZ4=0.0D0

      IF (DN(I).EQ.2) THEN
C         PRINT *,'PN1=2'
         M= (DVN(I)/DID(I))*4.0D0*COS(DPHI(I))*COS(DDELTA(I))
         DMX1= (DVN(I)/DID(I))*4.0D0*COS(DDELTA(I))*((V*DUX1 - U*DVX1)/V**2)
         DMX2= (DVN(I)/DID(I))*4.0D0*COS(DDELTA(I))*((V*DUX2 - U*DVX2)/V**2)
         DMX3= (DVN(I)/DID(I))*4.0D0*COS(DDELTA(I))*((V*DUX3 - U*DVX3)/V**2)
         DMX4= (DVN(I)/DID(I))*4.0D0*COS(DDELTA(I))*((V*DUX4 - U*DVX4)/V**2)
         DMY1= (DVN(I)/DID(I))*4.0D0*COS(DDELTA(I))*((V*DUY1 - U*DVY1)/V**2)
         DMY2= (DVN(I)/DID(I))*4.0D0*COS(DDELTA(I))*((V*DUY2 - U*DVY2)/V**2)
         DMY3= (DVN(I)/DID(I))*4.0D0*COS(DDELTA(I))*((V*DUY3 - U*DVY3)/V**2)
         DMY4= (DVN(I)/DID(I))*4.0D0*COS(DDELTA(I))*((V*DUY4 - U*DVY4)/V**2)
         DMZ1= (DVN(I)/DID(I))*4.0D0*COS(DDELTA(I))*((V*DUZ1 - U*DVZ1)/V**2)
         DMZ2= (DVN(I)/DID(I))*4.0D0*COS(DDELTA(I))*((V*DUZ2 - U*DVZ2)/V**2)
         DMZ3= (DVN(I)/DID(I))*4.0D0*COS(DDELTA(I))*((V*DUZ3 - U*DVZ3)/V**2)
         DMZ4= (DVN(I)/DID(I))*4.0D0*COS(DDELTA(I))*((V*DUZ4 - U*DVZ4)/V**2)
      ELSE IF (DN(I).EQ.3) THEN
C         PRINT *,'PN1=3'
         M= (DVN(I)/DID(I))*(9.0D0-12.0D0*(SIN(DPHI(I)))**2)*COS(DDELTA(I))
         DMX1= (DVN(I)/DID(I))*24.0D0*COS(DPHI(I))*COS(DDELTA(I))*((V*DUX1 - U*DVX1)/V**2)
         DMX2= (DVN(I)/DID(I))*24.0D0*COS(DPHI(I))*COS(DDELTA(I))*((V*DUX2 - U*DVX2)/V**2)
         DMX3= (DVN(I)/DID(I))*24.0D0*COS(DPHI(I))*COS(DDELTA(I))*((V*DUX3 - U*DVX3)/V**2)
         DMX4= (DVN(I)/DID(I))*24.0D0*COS(DPHI(I))*COS(DDELTA(I))*((V*DUX4 - U*DVX4)/V**2)
         DMY1= (DVN(I)/DID(I))*24.0D0*COS(DPHI(I))*COS(DDELTA(I))*((V*DUY1 - U*DVY1)/V**2)
         DMY2= (DVN(I)/DID(I))*24.0D0*COS(DPHI(I))*COS(DDELTA(I))*((V*DUY2 - U*DVY2)/V**2)
         DMY3= (DVN(I)/DID(I))*24.0D0*COS(DPHI(I))*COS(DDELTA(I))*((V*DUY3 - U*DVY3)/V**2)
         DMY4= (DVN(I)/DID(I))*24.0D0*COS(DPHI(I))*COS(DDELTA(I))*((V*DUY4 - U*DVY4)/V**2)
         DMZ1= (DVN(I)/DID(I))*24.0D0*COS(DPHI(I))*COS(DDELTA(I))*((V*DUZ1 - U*DVZ1)/V**2)
         DMZ2= (DVN(I)/DID(I))*24.0D0*COS(DPHI(I))*COS(DDELTA(I))*((V*DUZ2 - U*DVZ2)/V**2)
         DMZ3= (DVN(I)/DID(I))*24.0D0*COS(DPHI(I))*COS(DDELTA(I))*((V*DUZ3 - U*DVZ3)/V**2)
         DMZ4= (DVN(I)/DID(I))*24.0D0*COS(DPHI(I))*COS(DDELTA(I))*((V*DUZ4 - U*DVZ4)/V**2)
      ELSE IF (DN(I).EQ.4) THEN
C         PRINT *,'PN1=4'
         M= (DVN(I)/DID(I))*(32.0D0*(COS(DPHI(I)))**3 - 16.0D0*COS(DPHI(I)))*COS(DDELTA(I))
         DMX1= (DVN(I)/DID(I))*(96.0D0*(COS(DPHI(I))**2)-16.0D0)*COS(DDELTA(I))*((V*DUX1 - U*DVX1)/V**2)
         DMX2= (DVN(I)/DID(I))*(96.0D0*(COS(DPHI(I))**2)-16.0D0)*COS(DDELTA(I))*((V*DUX2 - U*DVX2)/V**2)
         DMX3= (DVN(I)/DID(I))*(96.0D0*(COS(DPHI(I))**2)-16.0D0)*COS(DDELTA(I))*((V*DUX3 - U*DVX3)/V**2)
         DMX4= (DVN(I)/DID(I))*(96.0D0*(COS(DPHI(I))**2)-16.0D0)*COS(DDELTA(I))*((V*DUX4 - U*DVX4)/V**2)
         DMY1= (DVN(I)/DID(I))*(96.0D0*(COS(DPHI(I))**2)-16.0D0)*COS(DDELTA(I))*((V*DUY1 - U*DVY1)/V**2)
         DMY2= (DVN(I)/DID(I))*(96.0D0*(COS(DPHI(I))**2)-16.0D0)*COS(DDELTA(I))*((V*DUY2 - U*DVY2)/V**2)
         DMY3= (DVN(I)/DID(I))*(96.0D0*(COS(DPHI(I))**2)-16.0D0)*COS(DDELTA(I))*((V*DUY3 - U*DVY3)/V**2)
         DMY4= (DVN(I)/DID(I))*(96.0D0*(COS(DPHI(I))**2)-16.0D0)*COS(DDELTA(I))*((V*DUY4 - U*DVY4)/V**2)
         DMZ1= (DVN(I)/DID(I))*(96.0D0*(COS(DPHI(I))**2)-16.0D0)*COS(DDELTA(I))*((V*DUZ1 - U*DVZ1)/V**2)
         DMZ2= (DVN(I)/DID(I))*(96.0D0*(COS(DPHI(I))**2)-16.0D0)*COS(DDELTA(I))*((V*DUZ2 - U*DVZ2)/V**2)
         DMZ3= (DVN(I)/DID(I))*(96.0D0*(COS(DPHI(I))**2)-16.0D0)*COS(DDELTA(I))*((V*DUZ3 - U*DVZ3)/V**2)
         DMZ4= (DVN(I)/DID(I))*(96.0D0*(COS(DPHI(I))**2)-16.0D0)*COS(DDELTA(I))*((V*DUZ4 - U*DVZ4)/V**2)
      END IF

      IF (DN2(I).EQ.1) THEN
C         PRINT *,'PN2=1'
         M= M+(DVN2(I)/DID(I))*COS(DDELTA2(I))
C
C DO NOTHING
C
      ELSE IF (DN2(I).EQ.2) THEN
C         PRINT *,'PN2=2'
         M= M+ (DVN2(I)/DID(I))*4.0D0*COS(DPHI(I))*COS(DDELTA(I))
         DMX1= DMX1+ (DVN2(I)/DID(I))*4.0D0*COS(DDELTA2(I))*((V*DUX1 - U*DVX1)/V**2)
         DMX2= DMX2+ (DVN2(I)/DID(I))*4.0D0*COS(DDELTA2(I))*((V*DUX2 - U*DVX2)/V**2)
         DMX3= DMX3+ (DVN2(I)/DID(I))*4.0D0*COS(DDELTA2(I))*((V*DUX3 - U*DVX3)/V**2)
         DMX4= DMX4+ (DVN2(I)/DID(I))*4.0D0*COS(DDELTA2(I))*((V*DUX4 - U*DVX4)/V**2)
         DMY1= DMY1+ (DVN2(I)/DID(I))*4.0D0*COS(DDELTA2(I))*((V*DUY1 - U*DVY1)/V**2)
         DMY2= DMY2+ (DVN2(I)/DID(I))*4.0D0*COS(DDELTA2(I))*((V*DUY2 - U*DVY2)/V**2)
         DMY3= DMY3+ (DVN2(I)/DID(I))*4.0D0*COS(DDELTA2(I))*((V*DUY3 - U*DVY3)/V**2)
         DMY4= DMY4+ (DVN2(I)/DID(I))*4.0D0*COS(DDELTA2(I))*((V*DUY4 - U*DVY4)/V**2)
         DMZ1= DMZ1+ (DVN2(I)/DID(I))*4.0D0*COS(DDELTA2(I))*((V*DUZ1 - U*DVZ1)/V**2)
         DMZ2= DMZ2+ (DVN2(I)/DID(I))*4.0D0*COS(DDELTA2(I))*((V*DUZ2 - U*DVZ2)/V**2)
         DMZ3= DMZ3+ (DVN2(I)/DID(I))*4.0D0*COS(DDELTA2(I))*((V*DUZ3 - U*DVZ3)/V**2)
         DMZ4= DMZ4+ (DVN2(I)/DID(I))*4.0D0*COS(DDELTA2(I))*((V*DUZ4 - U*DVZ4)/V**2)
      ELSE IF (DN2(I).EQ.3) THEN
C         PRINT *,'PN2=3'
         M= M+ (DVN2(I)/DID(I))*(9.0D0-12.0D0*(SIN(DPHI(I)))**2)*COS(DDELTA2(I))
         DMX1= DMX1+ (DVN2(I)/DID(I))*24.0D0*COS(DPHI(I))*COS(DDELTA2(I))*((V*DUX1 - U*DVX1)/V**2)
         DMX2= DMX2+ (DVN2(I)/DID(I))*24.0D0*COS(DPHI(I))*COS(DDELTA2(I))*((V*DUX2 - U*DVX2)/V**2)
         DMX3= DMX3+ (DVN2(I)/DID(I))*24.0D0*COS(DPHI(I))*COS(DDELTA2(I))*((V*DUX3 - U*DVX3)/V**2)
         DMX4= DMX4+ (DVN2(I)/DID(I))*24.0D0*COS(DPHI(I))*COS(DDELTA2(I))*((V*DUX4 - U*DVX4)/V**2)
         DMY1= DMY1+ (DVN2(I)/DID(I))*24.0D0*COS(DPHI(I))*COS(DDELTA2(I))*((V*DUY1 - U*DVY1)/V**2)
         DMY2= DMY2+ (DVN2(I)/DID(I))*24.0D0*COS(DPHI(I))*COS(DDELTA2(I))*((V*DUY2 - U*DVY2)/V**2)
         DMY3= DMY3+ (DVN2(I)/DID(I))*24.0D0*COS(DPHI(I))*COS(DDELTA2(I))*((V*DUY3 - U*DVY3)/V**2)
         DMY4= DMY4+ (DVN2(I)/DID(I))*24.0D0*COS(DPHI(I))*COS(DDELTA2(I))*((V*DUY4 - U*DVY4)/V**2)
         DMZ1= DMZ1+ (DVN2(I)/DID(I))*24.0D0*COS(DPHI(I))*COS(DDELTA2(I))*((V*DUZ1 - U*DVZ1)/V**2)
         DMZ2= DMZ2+ (DVN2(I)/DID(I))*24.0D0*COS(DPHI(I))*COS(DDELTA2(I))*((V*DUZ2 - U*DVZ2)/V**2)
         DMZ3= DMZ3+ (DVN2(I)/DID(I))*24.0D0*COS(DPHI(I))*COS(DDELTA2(I))*((V*DUZ3 - U*DVZ3)/V**2)
         DMZ4= DMZ4+ (DVN2(I)/DID(I))*24.0D0*COS(DPHI(I))*COS(DDELTA2(I))*((V*DUZ4 - U*DVZ4)/V**2)
      ELSE IF (DN2(I).EQ.4) THEN
C         PRINT *,'PN2=4'
         M= M+ (DVN2(I)/DID(I))*(32.0D0*(COS(DPHI(I)))**3 - 16.0D0*COS(DPHI(I)))*COS(DDELTA2(I))
         DMX1= DMX1+ (DVN2(I)/DID(I))*(96.0D0*(COS(DPHI(I))**2)-16.0D0)*COS(DDELTA2(I))*((V*DUX1 - U*DVX1)/V**2)
         DMX2= DMX2+ (DVN2(I)/DID(I))*(96.0D0*(COS(DPHI(I))**2)-16.0D0)*COS(DDELTA2(I))*((V*DUX2 - U*DVX2)/V**2)
         DMX3= DMX3+ (DVN2(I)/DID(I))*(96.0D0*(COS(DPHI(I))**2)-16.0D0)*COS(DDELTA2(I))*((V*DUX3 - U*DVX3)/V**2)
         DMX4= DMX4+ (DVN2(I)/DID(I))*(96.0D0*(COS(DPHI(I))**2)-16.0D0)*COS(DDELTA2(I))*((V*DUX4 - U*DVX4)/V**2)
         DMY1= DMY1+ (DVN2(I)/DID(I))*(96.0D0*(COS(DPHI(I))**2)-16.0D0)*COS(DDELTA2(I))*((V*DUY1 - U*DVY1)/V**2)
         DMY2= DMY2+ (DVN2(I)/DID(I))*(96.0D0*(COS(DPHI(I))**2)-16.0D0)*COS(DDELTA2(I))*((V*DUY2 - U*DVY2)/V**2)
         DMY3= DMY3+ (DVN2(I)/DID(I))*(96.0D0*(COS(DPHI(I))**2)-16.0D0)*COS(DDELTA2(I))*((V*DUY3 - U*DVY3)/V**2)
         DMY4= DMY4+ (DVN2(I)/DID(I))*(96.0D0*(COS(DPHI(I))**2)-16.0D0)*COS(DDELTA2(I))*((V*DUY4 - U*DVY4)/V**2)
         DMZ1= DMZ1+ (DVN2(I)/DID(I))*(96.0D0*(COS(DPHI(I))**2)-16.0D0)*COS(DDELTA2(I))*((V*DUZ1 - U*DVZ1)/V**2)
         DMZ2= DMZ2+ (DVN2(I)/DID(I))*(96.0D0*(COS(DPHI(I))**2)-16.0D0)*COS(DDELTA2(I))*((V*DUZ2 - U*DVZ2)/V**2)
         DMZ3= DMZ3+ (DVN2(I)/DID(I))*(96.0D0*(COS(DPHI(I))**2)-16.0D0)*COS(DDELTA2(I))*((V*DUZ3 - U*DVZ3)/V**2)
         DMZ4= DMZ4+ (DVN2(I)/DID(I))*(96.0D0*(COS(DPHI(I))**2)-16.0D0)*COS(DDELTA2(I))*((V*DUZ4 - U*DVZ4)/V**2)
      END IF

      IF (DN3(I).EQ.1) THEN
C         PRINT *,'PN3=1'
         M= M+(DVN3(I)/DID(I))*COS(DDELTA3(I))
C
C DO NOTHING
C
      ELSE IF (DN3(I).EQ.2) THEN
C         PRINT *,'PN3=2'
         M= M+ (DVN3(I)/DID(I))*4.0D0*COS(DPHI(I))*COS(DDELTA3(I))
         DMX1= DMX1+ (DVN3(I)/DID(I))*4.0D0*COS(DDELTA3(I))*((V*DUX1 - U*DVX1)/V**2)
         DMX2= DMX2+ (DVN3(I)/DID(I))*4.0D0*COS(DDELTA3(I))*((V*DUX2 - U*DVX2)/V**2)
         DMX3= DMX3+ (DVN3(I)/DID(I))*4.0D0*COS(DDELTA3(I))*((V*DUX3 - U*DVX3)/V**2)
         DMX4= DMX4+ (DVN3(I)/DID(I))*4.0D0*COS(DDELTA3(I))*((V*DUX4 - U*DVX4)/V**2)
         DMY1= DMY1+ (DVN3(I)/DID(I))*4.0D0*COS(DDELTA3(I))*((V*DUY1 - U*DVY1)/V**2)
         DMY2= DMY2+ (DVN3(I)/DID(I))*4.0D0*COS(DDELTA3(I))*((V*DUY2 - U*DVY2)/V**2)
         DMY3= DMY3+ (DVN3(I)/DID(I))*4.0D0*COS(DDELTA3(I))*((V*DUY3 - U*DVY3)/V**2)
         DMY4= DMY4+ (DVN3(I)/DID(I))*4.0D0*COS(DDELTA3(I))*((V*DUY4 - U*DVY4)/V**2)
         DMZ1= DMZ1+ (DVN3(I)/DID(I))*4.0D0*COS(DDELTA3(I))*((V*DUZ1 - U*DVZ1)/V**2)
         DMZ2= DMZ2+ (DVN3(I)/DID(I))*4.0D0*COS(DDELTA3(I))*((V*DUZ2 - U*DVZ2)/V**2)
         DMZ3= DMZ3+ (DVN3(I)/DID(I))*4.0D0*COS(DDELTA3(I))*((V*DUZ3 - U*DVZ3)/V**2)
         DMZ4= DMZ4+ (DVN3(I)/DID(I))*4.0D0*COS(DDELTA3(I))*((V*DUZ4 - U*DVZ4)/V**2)
      ELSE IF (DN3(I).EQ.3) THEN
C         PRINT *,'PN3=3'
         M= M+ (DVN3(I)/DID(I))*(9.0D0-12.0D0*(SIN(DPHI(I)))**2)*COS(DDELTA3(I))
         DMX1= DMX1+ (DVN3(I)/DID(I))*24.0D0*COS(DPHI(I))*COS(DDELTA3(I))*((V*DUX1 - U*DVX1)/V**2)
         DMX2= DMX2+ (DVN3(I)/DID(I))*24.0D0*COS(DPHI(I))*COS(DDELTA3(I))*((V*DUX2 - U*DVX2)/V**2)
         DMX3= DMX3+ (DVN3(I)/DID(I))*24.0D0*COS(DPHI(I))*COS(DDELTA3(I))*((V*DUX3 - U*DVX3)/V**2)
         DMX4= DMX4+ (DVN3(I)/DID(I))*24.0D0*COS(DPHI(I))*COS(DDELTA3(I))*((V*DUX4 - U*DVX4)/V**2)
         DMY1= DMY1+ (DVN3(I)/DID(I))*24.0D0*COS(DPHI(I))*COS(DDELTA3(I))*((V*DUY1 - U*DVY1)/V**2)
         DMY2= DMY2+ (DVN3(I)/DID(I))*24.0D0*COS(DPHI(I))*COS(DDELTA3(I))*((V*DUY2 - U*DVY2)/V**2)
         DMY3= DMY3+ (DVN3(I)/DID(I))*24.0D0*COS(DPHI(I))*COS(DDELTA3(I))*((V*DUY3 - U*DVY3)/V**2)
         DMY4= DMY4+ (DVN3(I)/DID(I))*24.0D0*COS(DPHI(I))*COS(DDELTA3(I))*((V*DUY4 - U*DVY4)/V**2)
         DMZ1= DMZ1+ (DVN3(I)/DID(I))*24.0D0*COS(DPHI(I))*COS(DDELTA3(I))*((V*DUZ1 - U*DVZ1)/V**2)
         DMZ2= DMZ2+ (DVN3(I)/DID(I))*24.0D0*COS(DPHI(I))*COS(DDELTA3(I))*((V*DUZ2 - U*DVZ2)/V**2)
         DMZ3= DMZ3+ (DVN3(I)/DID(I))*24.0D0*COS(DPHI(I))*COS(DDELTA3(I))*((V*DUZ3 - U*DVZ3)/V**2)
         DMZ4= DMZ4+ (DVN3(I)/DID(I))*24.0D0*COS(DPHI(I))*COS(DDELTA3(I))*((V*DUZ4 - U*DVZ4)/V**2)
      ELSE IF (DN3(I).EQ.4) THEN
C         PRINT *,'PN3=4'
         M= M+ (DVN3(I)/DID(I))*(32.0D0*(COS(DPHI(I)))**3 - 16.0D0*COS(DPHI(I)))*COS(DDELTA3(I))
         DMX1= DMX1+ (DVN3(I)/DID(I))*(96.0D0*(COS(DPHI(I))**2)-16.0D0)*COS(DDELTA3(I))*((V*DUX1 - U*DVX1)/V**2)
         DMX2= DMX2+ (DVN3(I)/DID(I))*(96.0D0*(COS(DPHI(I))**2)-16.0D0)*COS(DDELTA3(I))*((V*DUX2 - U*DVX2)/V**2)
         DMX3= DMX3+ (DVN3(I)/DID(I))*(96.0D0*(COS(DPHI(I))**2)-16.0D0)*COS(DDELTA3(I))*((V*DUX3 - U*DVX3)/V**2)
         DMX4= DMX4+ (DVN3(I)/DID(I))*(96.0D0*(COS(DPHI(I))**2)-16.0D0)*COS(DDELTA3(I))*((V*DUX4 - U*DVX4)/V**2)
         DMY1= DMY1+ (DVN3(I)/DID(I))*(96.0D0*(COS(DPHI(I))**2)-16.0D0)*COS(DDELTA3(I))*((V*DUY1 - U*DVY1)/V**2)
         DMY2= DMY2+ (DVN3(I)/DID(I))*(96.0D0*(COS(DPHI(I))**2)-16.0D0)*COS(DDELTA3(I))*((V*DUY2 - U*DVY2)/V**2)
         DMY3= DMY3+ (DVN3(I)/DID(I))*(96.0D0*(COS(DPHI(I))**2)-16.0D0)*COS(DDELTA3(I))*((V*DUY3 - U*DVY3)/V**2)
         DMY4= DMY4+ (DVN3(I)/DID(I))*(96.0D0*(COS(DPHI(I))**2)-16.0D0)*COS(DDELTA3(I))*((V*DUY4 - U*DVY4)/V**2)
         DMZ1= DMZ1+ (DVN3(I)/DID(I))*(96.0D0*(COS(DPHI(I))**2)-16.0D0)*COS(DDELTA3(I))*((V*DUZ1 - U*DVZ1)/V**2)
         DMZ2= DMZ2+ (DVN3(I)/DID(I))*(96.0D0*(COS(DPHI(I))**2)-16.0D0)*COS(DDELTA3(I))*((V*DUZ2 - U*DVZ2)/V**2)
         DMZ3= DMZ3+ (DVN3(I)/DID(I))*(96.0D0*(COS(DPHI(I))**2)-16.0D0)*COS(DDELTA3(I))*((V*DUZ3 - U*DVZ3)/V**2)
         DMZ4= DMZ4+ (DVN3(I)/DID(I))*(96.0D0*(COS(DPHI(I))**2)-16.0D0)*COS(DDELTA3(I))*((V*DUZ4 - U*DVZ4)/V**2)
      END IF


      RETURN

      END


      SUBROUTINE HAIRYIMPHELL(U,V)
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE
      DOUBLE PRECISION      DUX1,DUX2,DUX3,DUX4,DUY1,DUY2,DUY3,DUY4,DUZ1,DUZ2,DUZ3,DUZ4
      COMMON /DU/           DUX1,DUX2,DUX3,DUX4,DUY1,DUY2,DUY3,DUY4,DUZ1,DUZ2,DUZ3,DUZ4
      DOUBLE PRECISION      DVX1,DVX2,DVX3,DVX4,DVY1,DVY2,DVY3,DVY4,DVZ1,DVZ2,DVZ3,DVZ4
      COMMON /DV/           DVX1,DVX2,DVX3,DVX4,DVY1,DVY2,DVY3,DVY4,DVZ1,DVZ2,DVZ3,DVZ4
      DOUBLE PRECISION      DMX1,DMX2,DMX3,DMX4,DMY1,DMY2,DMY3,DMY4,DMZ1,DMZ2,DMZ3,DMZ4
      COMMON /DM/           DMX1,DMX2,DMX3,DMX4,DMY1,DMY2,DMY3,DMY4,DMZ1,DMZ2,DMZ3,DMZ4
      DOUBLE PRECISION      U,V
      DOUBLE PRECISION      M,N
      COMMON /MN/           M,N

C      PRINT *,'IN HAIRYIMPHELL',I
      M=0.0D0
      DMX1=0.0D0
      DMX2=0.0D0
      DMX3=0.0D0
      DMX4=0.0D0
      DMY1=0.0D0
      DMY2=0.0D0
      DMY3=0.0D0
      DMY4=0.0D0
      DMZ1=0.0D0
      DMZ2=0.0D0
      DMZ3=0.0D0
      DMZ4=0.0D0

      IF (IN1(I).EQ.2) THEN
C         PRINT *,'PN1=2'
         M= IVN(I)*4.0D0*COS(IPHI(I))*COS(IDELTA(I))
         DMX1= IVN(I)*4.0D0*COS(IDELTA(I))*((V*DUX1 - U*DVX1)/V**2)
         DMX2= IVN(I)*4.0D0*COS(IDELTA(I))*((V*DUX2 - U*DVX2)/V**2)
         DMX3= IVN(I)*4.0D0*COS(IDELTA(I))*((V*DUX3 - U*DVX3)/V**2)
         DMX4= IVN(I)*4.0D0*COS(IDELTA(I))*((V*DUX4 - U*DVX4)/V**2)
         DMY1= IVN(I)*4.0D0*COS(IDELTA(I))*((V*DUY1 - U*DVY1)/V**2)
         DMY2= IVN(I)*4.0D0*COS(IDELTA(I))*((V*DUY2 - U*DVY2)/V**2)
         DMY3= IVN(I)*4.0D0*COS(IDELTA(I))*((V*DUY3 - U*DVY3)/V**2)
         DMY4= IVN(I)*4.0D0*COS(IDELTA(I))*((V*DUY4 - U*DVY4)/V**2)
         DMZ1= IVN(I)*4.0D0*COS(IDELTA(I))*((V*DUZ1 - U*DVZ1)/V**2)
         DMZ2= IVN(I)*4.0D0*COS(IDELTA(I))*((V*DUZ2 - U*DVZ2)/V**2)
         DMZ3= IVN(I)*4.0D0*COS(IDELTA(I))*((V*DUZ3 - U*DVZ3)/V**2)
         DMZ4= IVN(I)*4.0D0*COS(IDELTA(I))*((V*DUZ4 - U*DVZ4)/V**2)
      ELSE 
C         PRINT *,'EEEEK! COCK UP IN HAIRYIMPHELL!!!'
      END IF


      RETURN

      END

      SUBROUTINE AMBERDUMP(ACOORDS,AFILENAME)
      USE COMMONS
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE
      DOUBLE PRECISION  ACOORDS(3*NATOMS)
      CHARACTER(LEN=*)      AFILENAME
      INTEGER           DUMPINT

      OPEN(UNIT=19,FILE=AFILENAME,STATUS='UNKNOWN')
      DO DUMPINT=1,ATOMS
         WRITE (19,FMT='(A1,2X,A2,2X,I3,2X,I3,2X,F20.10,X,F20.10,X,F20.10)')
     1      LABEL(DUMPINT),TYPECH(DUMPINT),DUMPINT,BONDEDTO(DUMPINT),
     2      ACOORDS(3*DUMPINT-2),ACOORDS(3*DUMPINT-1),ACOORDS(3*DUMPINT)
      END DO
      WRITE (19,FMT='(A)') 'END'
      WRITE (19,FMT='(A)') ' '
      WRITE (19,FMT='(A4,7X,I2)') 'LOOP',RINGS
      DO DUMPINT=1,RINGS
         WRITE (19,FMT='(I3,2X,I3)') LOOPATOM(2*DUMPINT-1),LOOPATOM(2*DUMPINT)
      END DO
      WRITE (19,FMT='(A)') ' '
      WRITE (19,FMT='(A)') 'CHARGES'
      DO DUMPINT=1,ATOMS
         WRITE (19,FMT='(I3,2X,F7.4)') DUMPINT,PQ(DUMPINT)/18.2223
      END DO
      WRITE (19,FMT='(A)') 'END'
      CLOSE (19)


      RETURN

      END


      SUBROUTINE ACONNECTDUMP
      USE MODAMBER
      USE MODAMBER2
      IMPLICIT NONE
      INTEGER       CI,CJ

      OPEN (UNIT=17,FILE='AMBER_CONNECTIVITY',STATUS='UNKNOWN')
      WRITE (17,FMT='(A)') 'ANGLES'
      DO CI=1,ANG
         WRITE (17,FMT='(I3,A,I3,A,I3)') AA1(CI),'-',AA2(CI),'-',AA3(CI)
      END DO
      WRITE (17,FMT='(A)') 'TORSIONS'
      DO CI=1,T
         WRITE (17,FMT='(I3,A,I3,A,I3,A,I3,A,4F10.5)') DA1(CI),'-',DA2(CI),'-',DA3(CI),'-',DA4(CI),'....',DVN(CI),DID(CI),DDELTA(CI)
     1         ,DN(CI)
      END DO
      WRITE (17,FMT='(A)') 'IMPROPERS'
      DO CI=1,IMP
         WRITE (17,FMT='(I3,A,I3,A,I3,A,I3)') IA1(CI),'-',IA2(CI),'-',IA3(CI),'-',IA4(CI)
      END DO
      WRITE (17,FMT='(A)') 'PAIRWISE PARAMETERS'
      DO CI=1,ATOMS
        DO CJ=1,ATOMS
          WRITE (17,FMT='(I3,2X,I3,2X,I1,2X,I1,2X,I1,2X,F20.7,2X,F20.7)') 
     1          CI,CJ,BONDS(CI,CJ),
     2          ONE_THREE(CI,CJ),ONE_FOUR(CI,CJ),VDWA(CI,CJ),VDWB(CI,CJ)
        END DO
      END DO
      CLOSE (17)


      RETURN

      END
