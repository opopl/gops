
      SUBROUTINE POTENTIAL(X,GRAD,EREAL,GRADT,SECT)
!op226> Declarations {{{ 
      USE COMMONS
      USE PORFUNCS
      USE BLN

      IMPLICIT NONE
      
      LOGICAL GRADT, FTEST, SECT, EVAP, COMPON, YESNO, GUIDET, evapreject
      COMMON /CO/ COMPON
      INTEGER J1, J2, J3, NPCALL, PERM(NATOMS), NPERM, NORBIT1, NORBIT2, CSMIT
      CHARACTER FNAME*80, DUMM*4
      DOUBLE PRECISION EREAL, GRAD(*), X(*), DUMMY2, GEMAX, XG, YG, ZG, RMAT(3,3), XD, YD, ZD, XTEMP(3*NATOMS),
     1                 GRADLJ(3*NATOMS), EREALLJ, GRADMF(3*NATOMS), EREALMF, TERMLJ, TERMMF, GRADDUM1(3*NATOMS), AVVAL,
     2                 GRADDUM2(3*NATOMS), EDUM1, EDUM2, DUMMY(3*NATOMS), DIST2, WORSTRAD, GBDUMMY, QE, QX, DISTANCE, AA(3),
     3                 SAVECSMNORM, CSMRMS, CSMGRAD(3), SAVECSMPMAT(3,3), SAVECSMIMAGES(3*NATOMS*CSMGPINDEX)
C    2                 GRADNUM(3*NATOMS), DISP, V1, V2
C     DOUBLE PRECISION h0(NATOMS*3,NATOMS*3), h1(NATOMS*3,NATOMS*3), 
C    1                 ee(NATOMS*3), ev(NATOMS*3,NATOMS*3), w(NATOMS*3,NATOMS*3) 
      CHARACTER(LEN=3) CSMGPSAVE
      INTEGER CSMGPINDEXSAVE
      DOUBLE PRECISION PTGPSAVE(3,3,2*CSMGPINDEX), CSMNORMSAVE

      LOGICAL SOCOUPLE, GUIDECHANGET, CSMDOGUIDET
      COMMON /FAIL/ FTEST
      COMMON /EV/ EVAP, EVAPREJECT
      COMMON /GD/ GUIDECHANGET, GUIDET, CSMDOGUIDET
      COMMON /PCALL/ NPCALL
      COMMON /CSMAVVAL/ AVVAL, CSMRMS, CSMIT
      INTEGER BRUN

      INTEGER NQTOT
      COMMON /TOT/ NQTOT

      DOUBLE PRECISION EPLUS, EMINUS, GRADDUM(3*NATOMS), DIFF

!|gd351>
!only needed for numerical derivative checks
!      DOUBLE PRECISION :: TEMPX, TEMPLEFT, TEMPRIGHT
!      DOUBLE PRECISION :: NUMGRAD(3*NATOMS), ATEMP(6)
!      INTEGER :: I, I0, COUNTGOOD, COUNTBAD, ITEMP(6)

!<gd351|

Cop226> }}}
!op226>{{{ 
      GUIDECHANGET=.FALSE.
     
10    CONTINUE

      ELSE IF (P46) THEN
         CALL P46MERDIFF(X,NATOMS,GRAD,EREAL,GRADT)
      ELSE IF (G46) THEN
         CALL G46MERDIFF(X,NATOMS,GRAD,EREAL,GRADT)
      ELSE IF (BLNT) THEN
         CALL BLN(X,GRAD,EREAL,GRADT)
      ELSE
         CALL RAD(X,GRAD,EREAL,GRADT)
         IF (EVAPREJECT) return
         IF (CUTT) THEN
            CALL LJCUT(X,GRAD,EREAL,GRADT,SECT)
         ELSE
            CALL LJ(X,GRAD,EREAL,GRADT,SECT)
         ENDIF
      ENDIF

C  --------------- End of possible potentials - now add fields if required ------------------------------

      IF (PULLT) THEN
         EREAL=EREAL-PFORCE*(X(3*(PATOM1-1)+3)-X(3*(PATOM2-1)+3))
         GRAD(3*(PATOM1-1)+3)=GRAD(3*(PATOM1-1)+3)-PFORCE
         GRAD(3*(PATOM2-1)+3)=GRAD(3*(PATOM2-1)+3)+PFORCE
      ENDIF
      RETURN
      END
