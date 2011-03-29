! ALL THIS IS CONCERNING THE ALIGNMENT OF START AND FINISH FOR AMBER (ALL ATOM)
! AIM IS TO FIX "CHIRALITY" AROUND EVERY X - CH_2 - Y CARBON SO THAT 
! INTERPOLATION IS POSSIBLE
! MSB50

      SUBROUTINE FINDCHIRALH(PTEST)
      USE COMMONS,ONLY: NATOMS
      USE MODAMBER9
      USE INTCOMMONS, ONLY: MBBONDNUM, MBADJACENT 

      ! PROCHIRALH - IS THIS ATOM NEXT TO A PROCHIRAL ONE - LEN NATOMS, BUT CAN 
      ! BE TRUE ONLY FOR HYDROGENS
      ! PROCHIRALCNT - CENTREATOM WHICH IS PROCHIRAL
      ! PROCHIRALCNT(H1) = PROCHIRALCNT(H2) = C
      ! PROCHIRALCNT(C) - 4TH ATOM TO CALCULATE IMPROPER TO FIX PROCHIRALITY
      ! PROCHIRALCNT(PROCHIRAL_CENTRE) = 4TH ATOM (X OR Y, WHICHEVER HAS
      !                                 LOWEST ATOM INDEX)
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: PTEST
      INTEGER :: ATM1, ATM2, RNUM,I, CNTR, B, A
      CHARACTER(LEN=4):: PHERING(6), PRORING(5), TRPRING5(5),  &
     &     TRPRING6(6), TRPFRNG(6), HISRING(5), RES_RNUM
      LOGICAL :: RINGATM(NATOMS)
      INTEGER :: PROCH(2), HYD, OTHER, NONHYD(4)

      NBOND = NBONA + NBONH
      IF(.NOT.ALLOCATED(IBIB)) THEN
         ALLOCATE(IBIB(NBONA+NBONH))
         DO I =1,NBONA
            IBIB(I) = IX(IIBA+I-1)/3+1
         ENDDO
         DO I=1, NBONH
            IBIB(NBONA+I) = IX(IIBH+I-1)/3+1
         ENDDO
      ENDIF

      IF(.NOT.ALLOCATED(JBJB)) THEN
          ALLOCATE(JBJB(NBONA+NBONH))
          DO I =1,NBONA
             JBJB(I) = IX(IJBA+I-1)/3+1
          ENDDO
          DO I=1, NBONH
             JBJB(NBONA+I) = IX(IJBH+I-1)/3+1
          ENDDO
      ENDIF

      IF (.NOT.ALLOCATED(MBBONDNUM)) ALLOCATE(MBBONDNUM(NATOMS))
      IF (.NOT.ALLOCATED(MBADJACENT)) ALLOCATE(MBADJACENT(NATOMS,6))
      MBBONDNUM(:)=0

      ! RING DEFINITIONS FOR KNOWN RING RESIDUES
      PHERING = (/'CG  ', 'CD1 ', 'CE1 ', 'CZ  ', 'CE2 ', 'CD2 '/)
      PRORING = (/'N   ', 'CA  ', 'CB  ', 'CG  ', 'CD  '/)
      TRPRING5 = (/'CG  ', 'CD2 ', 'CE2 ', 'NE1 ', 'CD1 '/)
      TRPRING6 = (/'CD2 ', 'CE3 ', 'CZ3 ', 'CH2 ', 'CZ2 ', 'CE2 '/)
      TRPFRNG = (/'CE2 ', 'CD2 ', 'NE1 ', 'CG  ', 'CZ2 ', 'CE3 '/)
      HISRING = (/'CG  ', 'ND1 ', 'CE1 ', 'NE2 ', 'CD2 '/)

      DO B=1,NBOND
         ATM1 = IBIB(B); ATM2 = JBJB(B) 
         MBBONDNUM(ATM1) = MBBONDNUM(ATM1) + 1
         MBBONDNUM(ATM2) = MBBONDNUM(ATM2) + 1

         IF (MBBONDNUM(ATM1).GT.6.OR.MBBONDNUM(ATM2).GT.6) THEN
            PRINT*, 'GETNATINT>> ERROR! ATOM HAS MORE THAN 6 NEIGHBORS!'
            PRINT*, ATM1, MBBONDNUM(ATM1), ATM2, MBBONDNUM(ATM2)
            STOP
         ENDIF

         MBADJACENT(ATM1,MBBONDNUM(ATM1)) = ATM2
         MBADJACENT(ATM2,MBBONDNUM(ATM2)) = ATM1
      ENDDO

      RINGATM(:) = .FALSE.
      DO RNUM = 1, NRES
         ! PULL OUT RINGS FROM THE RING RESIDUES
         ! LABEL RING ATOMS IN RINGATM LIST
         RES_RNUM = IH(M02+RNUM-1)
         IF (RES_RNUM.EQ.'PHE'.OR.RES_RNUM.EQ.'TYR' &
     &      .OR.RES_RNUM.EQ.'NPHE'.OR.RES_RNUM.EQ.'NTYR' &
     &      .OR.RES_RNUM.EQ.'CPHE'.OR.RES_RNUM.EQ.'CTYR') THEN
            DO A = 1,6
               CALL AMB_PATOM(ATM1, RNUM, PHERING(A))
               IF (ATM1.LT.0) THEN
                  PRINT*, 'ERROR! COULD NOT FIND ATOM: ', PHERING(A), RNUM
                  STOP
               ENDIF
               RINGATM(ATM1) = .TRUE.
            ENDDO
         ELSE IF (RES_RNUM.EQ.'HIE'.OR.RES_RNUM.EQ.'HID'.OR.  &
     &          RES_RNUM.EQ.'HSP'    &
     &          .OR.RES_RNUM.EQ.'NHIE'.OR.RES_RNUM.EQ.'NHID'.OR. &
     &          RES_RNUM.EQ.'NHSP'.OR.RES_RNUM.EQ.'CHIE'.OR.RES_RNUM &
     &          .EQ.'CHID'.OR. RES_RNUM.EQ.'CHSP') THEN
            DO A = 1,5
               CALL AMB_PATOM(ATM1, RNUM, HISRING(A))
               IF (ATM1.LT.0) THEN
                  PRINT*, 'ERROR! COULD NOT FIND ATOM: ', HISRING(A), RNUM
                  STOP
               ENDIF
               RINGATM(ATM1) = .TRUE.
            ENDDO
         ELSE IF ((RES_RNUM.EQ.'PRO'.OR.RES_RNUM.EQ.'NPRO' &
     &        .OR. RES_RNUM.EQ.'CPRO')) THEN
            DO A = 1,5
               CALL AMB_PATOM(ATM1, RNUM, PRORING(A))
               IF (ATM1.LT.0) THEN
                  PRINT*, 'ERROR! COULD NOT FIND ATOM: ', PRORING(A), RNUM
                  STOP
               ENDIF
               RINGATM(ATM1) = .TRUE.
            ENDDO
         ELSE IF (RES_RNUM.EQ.'TRP'.OR.RES_RNUM.EQ.'NTRP'&
     &      .OR.RES_RNUM.EQ.'CTRP') THEN
            DO A = 1,5
               CALL AMB_PATOM(ATM1, RNUM, TRPRING5(A))
               IF (ATM1.LT.0) THEN
                  PRINT*, 'ERROR! COULD NOT FIND ATOM: ', TRPRING5(A), RNUM
                  STOP
               ENDIF
               RINGATM(ATM1) = .TRUE.
            ENDDO
            DO A = 1,6
               CALL AMB_PATOM(ATM1, RNUM, TRPRING6(A))
               IF (ATM1.LT.0) THEN
                  PRINT*, 'ERROR! COULD NOT FIND ATOM: ', TRPRING5(A), RNUM
                  STOP
               ENDIF
               RINGATM(ATM1) = .TRUE.
            ENDDO
        ENDIF
      ENDDO

      !FIND PROCHIRAL CENTRES: 4 BONDS, EXACTLY 2 H'S 
      IF (.NOT.ALLOCATED(PROCHIRALH)) ALLOCATE(PROCHIRALH(NATOMS))
      IF (.NOT.ALLOCATED(PROCHIRALCNT)) ALLOCATE(PROCHIRALCNT(NATOMS))
      PROCHIRALH(:)=.FALSE.
      DO CNTR=1, NATOMS
         IF (RINGATM(CNTR)) CYCLE
         IF (MBBONDNUM(CNTR).EQ.4) THEN
            HYD = 0; OTHER = 0
            DO I = 1, MBBONDNUM(CNTR)
               ATM1=MBADJACENT(CNTR,I)
               IF (IH(M04+ATM1-1)(:1).EQ."H") THEN
                  HYD = HYD + 1
                  PROCH(HYD)=ATM1
               ELSE
                  OTHER = OTHER + 1
                  NONHYD(OTHER) = ATM1
               ENDIF
            ENDDO       
            IF (HYD == 2) THEN
               PROCHIRALH(PROCH(1)) = .TRUE.
               PROCHIRALH(PROCH(2)) = .TRUE.
               PROCHIRALCNT(PROCH(1)) = CNTR
               PROCHIRALCNT(PROCH(2)) = CNTR
               IF (NONHYD(1) .LT. NONHYD(2)) THEN
                 PROCHIRALCNT(CNTR) = NONHYD(1)
               ELSE
                 PROCHIRALCNT(CNTR) = NONHYD(2) 
                ENDIF 
            ENDIF
         ENDIF
      ENDDO

        !DO I=1,NATOMS
        !  PRINT*, I, "PROCHIRAL", PROCHIRALH(I), PROCHIRALCNT(I)
        !ENDDO

      END SUBROUTINE FINDCHIRALH 


!***********************************************************************
      SUBROUTINE ALLCHIRALH_ALIGN(RS,RF,ATM1,ATM2,PTEST)
      USE COMMONS, ONLY: NATOMS
      USE MODAMBER9      
      USE INTCOMMONS, ONLY : INTMINPERMT
      IMPLICIT NONE
      
      DOUBLE PRECISION, INTENT(IN)    :: RS(3*NATOMS)
      DOUBLE PRECISION, INTENT(INOUT) :: RF(3*NATOMS)
      LOGICAL, INTENT(IN)             :: PTEST
      ! ATOMS WHICH HAVE TO BE PERMUTED MAINTAINING PROCHIRALITY 
      INTEGER, INTENT(IN)             :: ATM1, ATM2
      DOUBLE PRECISION                :: DIHED1, DIHED2
      DOUBLE PRECISION                :: TMP(3)
      INTEGER                         :: I1, J1
      DOUBLE PRECISION                :: DISTANCE1, DISTANCE2, INTDIST
      DOUBLE PRECISION                :: COORDS_COPY(3*NATOMS), RSTART(3*NATOMS)
      DOUBLE PRECISION                ::  TMP2(3), TMP3(3)
 
      I1 = PROCHIRALCNT(ATM1)
      J1 = PROCHIRALCNT(I1)
      CALL AMBERDIHEDR(RS, NATOMS,ATM1,ATM2,I1,J1, DIHED1)
      CALL AMBERDIHEDR(RF, NATOMS,ATM1,ATM2,I1,J1, DIHED2)
      IF (PTEST) PRINT*, "ALLCHIRALH, DIHED1, DIHED2", DIHED1, DIHED2
      IF (DIHED1*DIHED2.LT.0) THEN !SWAP RF HYDROGENS
         IF (PTEST) PRINT*, "ALLCHIRALH> SWAPP", ATM1, ATM2
         TMP(:) = RF(3*(ATM1-1)+1:3*ATM1)
         RF(3*(ATM1-1)+1:3*ATM1) = RF(3*(ATM2-1)+1:3*ATM2)
         RF(3*(ATM2-1)+1:3*ATM2) = TMP(:)
      ENDIF

 
      END SUBROUTINE ALLCHIRALH_ALIGN

!*************************************************************************   

! ***********************************************************************
      SUBROUTINE GLYINTPERM(RS, RF, PTEST)
      USE INTCOMMONS, ONLY: NDIH, PERMNEIGHBOURS, PERMCHAIN, NGLYDIH, GLYDIH !MSB50
      USE COMMONS, ONLY: NATOMS
      IMPLICIT NONE
    
      DOUBLE PRECISION, INTENT(IN)    :: RS(3*NATOMS)
      DOUBLE PRECISION, INTENT(INOUT) :: RF(3*NATOMS)
      LOGICAL, INTENT(IN)             :: PTEST
       
      INTEGER          :: III, JJJ
      DOUBLE PRECISION :: GLYSTART(4), GLYFIN(4), GLYDIFF(4), GLYDIST1, GLYDIST
      DOUBLE PRECISION :: DIHED1, RFTMP(3*NATOMS)
      INTEGER          :: I1, I2, I3, I4, PLUS

      RFTMP(:) = RF(:)
      DO III=1,NGLYDIH
         DO JJJ =2,3
            I1=GLYDIH(III,JJJ+3)
            I2=GLYDIH(III,JJJ)
            I3=GLYDIH(III,1)
            CALL AMBERDIHEDR(RS,NATOMS,I1,I2,I3,GLYDIH(III,7), DIHED1)
            GLYSTART(2*(JJJ-1)-1)=DIHED1
            CALL AMBERDIHEDR(RS,NATOMS,I1,I2,I3,GLYDIH(III,8), DIHED1)
            GLYSTART(2*(JJJ-1))=DIHED1
            CALL AMBERDIHEDR(RF,NATOMS,I1,I2,I3,GLYDIH(III,7), GLYFIN(2*(JJJ-1)-1))
            CALL AMBERDIHEDR(RF,NATOMS,I1,I2,I3,GLYDIH(III,8), GLYFIN(2*(JJJ-1)))
         ENDDO
         GLYDIFF(:)=GLYFIN(:)-GLYSTART(:)
         DO JJJ=1,4
            IF (GLYDIFF(JJJ).GT.180.0D0) THEN
                GLYDIFF(JJJ)=360.0D0-GLYDIFF(JJJ)
            ELSEIF (GLYDIFF(JJJ).LT.-180.0D0) THEN
                GLYDIFF(JJJ)=360+GLYDIFF(JJJ)
            ENDIF
         ENDDO
         GLYDIST1 = DOT_PRODUCT(GLYDIFF,GLYDIFF)
         IF (PTEST) PRINT*, "GLYDIST", GLYDIST1
         RF(3*(GLYDIH(III,7)-1)+1:3*GLYDIH(III,7)) = RFTMP(3*(GLYDIH(III,8)-1)+1:3*GLYDIH(III,8))
         RF(3*(GLYDIH(III,8)-1)+1:3*GLYDIH(III,8)) = RFTMP(3*(GLYDIH(III,7)-1)+1:3*GLYDIH(III,7))
         DO JJJ =2,3
            I1=GLYDIH(III,JJJ+3)
            I2=GLYDIH(III,JJJ)
            I3=GLYDIH(III,1)
            CALL AMBERDIHEDR(RF,NATOMS,I1,I2,I3,GLYDIH(III,7), GLYFIN(2*(JJJ-1)-1))
            CALL AMBERDIHEDR(RF,NATOMS,I1,I2,I3,GLYDIH(III,8), GLYFIN(2*(JJJ-1)))
         ENDDO
         GLYDIFF(:)=GLYFIN(:)-GLYSTART(:)
         DO JJJ=1,4
            IF (GLYDIFF(JJJ).GT.180.0D0) THEN
                GLYDIFF(JJJ)=360.0D0-GLYDIFF(JJJ)
            ELSEIF (GLYDIFF(JJJ).LT.-180.0D0) THEN
                GLYDIFF(JJJ)=360+GLYDIFF(JJJ)
            ENDIF
         ENDDO
         GLYDIST = DOT_PRODUCT(GLYDIFF,GLYDIFF)
         IF (PTEST) PRINT*, "GLYDIST", GLYDIST
         IF (GLYDIST.LT.GLYDIST1) THEN
            IF (PTEST) PRINT*, "KEEP PERMUTATION"
            !RFMP(:) = RF(:)
         ELSE
            RF(:)=RFTMP(:)
         ENDIF
      ENDDO

      END SUBROUTINE GLYINTPERM

! **************************************************************************
      SUBROUTINE MINPERM_CHIRAL(RS, RF,DISTANCE,RMAT,PTEST)
      !MINIMISE CARTESIAN DISTANCE WITHOUT PERMUTING THEN
      !PERMUTE AND CALCULATE CART DISTANCE USING A DUMB DOT PRODUCT
      !MINPERM WOULD BE BETTER, BUT HERE PROCHIRAL METHYL CARBONS ARE
      !FIXED TO HAVE SAME CHIRALITY
      USE KEY, ONLY : NPERMGROUP, NPERMSIZE, PERMGROUP, NSETS, SETS, NABT, AMBERT
      USE MODAMBER9, ONLY:NRES, IH, M02,  IX, I02, M04 , PROCHIRALH
      USE COMMONS, ONLY : NATOMS
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN) :: RS(3*NATOMS)
      DOUBLE PRECISION, INTENT(INOUT) :: RF(3*NATOMS), RMAT(3,3)
      DOUBLE PRECISION, INTENT(OUT) :: DISTANCE
      LOGICAL, INTENT(IN) :: PTEST

      INTEGER :: PG, SW, I, J, K, A1, A2, A3, B1, B2, B3, START,III
      DOUBLE PRECISION :: CURDIST, DIST 
      DOUBLE PRECISION :: RFTMP(3*NATOMS)
      DOUBLE PRECISION :: TMP(3), TMP2(3), TMP3(3)

      CURDIST = SQRT(DOT_PRODUCT(RF-RS, RF-RS))
      START = 1
      DO PG = 1,NPERMGROUP

        IF (NPERMSIZE(PG).EQ.2) THEN
           RFTMP(:) = RF(:)
           ! A1 AND A2 ARE THE ATOM NUMBERS TO PERMUTE
           A1 = PERMGROUP(START)
           A2 = PERMGROUP(START+1) 
 
           IF (.NOT.PROCHIRALH(A1)) THEN 

              IF (PTEST) PRINT*, 'SWAPPING ATOMS: ', A1, A2
              RF(3*(A1-1)+1:3*A1) = RFTMP(3*(A2-1)+1:3*A2)
              RF(3*(A2-1)+1:3*A2) = RFTMP(3*(A1-1)+1:3*A1)
            
              !MOVE ANY OTHER GROUPS THAT MUST BE DRAGGED ALONG
              DO SW = 1,NSETS(PG)
!                A1 = SWAP1(PG,SW)
!                A2 = SWAP2(PG,SW)
                 A1= SETS(PERMGROUP(START),SW)
                 A2= SETS(PERMGROUP(START+1),SW)
            
                 IF (PTEST) PRINT*, 'DRAGGING SWAP: ', A1, A2
            
                 RF(3*(A1-1)+1:3*A1) = RFTMP(3*(A2-1)+1:3*A2)
                 RF(3*(A2-1)+1:3*A2) = RFTMP(3*(A1-1)+1:3*A1)
              ENDDO
           ELSE
               CALL ALLCHIRALH_ALIGN(RS,RF,A1,A2,PTEST)
           ENDIF

           DIST = SQRT(DOT_PRODUCT(RF-RS, RF-RS))
           IF (PTEST) PRINT*, "DISTANCES", DIST, "CURRENT", CURDIST

           IF (DIST.GT.CURDIST.AND..NOT.PROCHIRALH(A1)) THEN
              IF (PTEST) PRINT*, "UNDO PERMUTATION"
              RF(:) = RFTMP(:)
           ELSEIF (DIST.LT.CURDIST) THEN
              IF (PTEST) PRINT*, "KEEP PERMUTATION"
              CURDIST = DIST 
           ELSEIF (PROCHIRALH(A1)) THEN
              CURDIST = DIST
           ENDIF

        ELSE IF (NPERMSIZE(PG).EQ.3) THEN
          A1 = PERMGROUP(START)
          A2 = PERMGROUP(START+1)
          A3 = PERMGROUP(START+2)

          TMP(:) = RF(3*(A1-1)+1:3*A1)
          TMP2(:) = RF(3*(A2-1)+1:3*A2)
          TMP3(:) = RF(3*(A3-1)+1:3*A3)

          RFTMP(:) = RF(:)
          DO I = 0,2
             DO J = 0,2
                IF (J.EQ.I.OR.(I.EQ.0.AND.J.EQ.1)) CYCLE
                DO K = 0,2
                   IF (K.EQ.I.OR.K.EQ.J) CYCLE

                   B1 = PERMGROUP(START+I)
                   B2 = PERMGROUP(START+J)
                   B3 = PERMGROUP(START+K)

                   IF (PTEST) PRINT*, 'PERMUTING: ', B1, B2, B3

                   RF(3*(B1-1)+1:3*B1) = TMP(:)
                   RF(3*(B2-1)+1:3*B2) = TMP2(:)
                   RF(3*(B3-1)+1:3*B3) = TMP3(:)

                   DIST = SQRT(DOT_PRODUCT(RF-RS, RF-RS))

                   IF (PTEST) PRINT*, 'NEW DISTANCE: ', DIST

                   IF (DIST.LT.CURDIST) THEN
                      RFTMP(:) = RF(:)
                      CURDIST = DIST
                      IF (PTEST) PRINT*, 'KEEP PERMUTATION'
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
          RF(:) = RFTMP(:)
       ELSE
          PRINT*, 'ERROR! INTMINPERM IS ONLY SET UP FOR PERMUTATION GROUPS &
               &          OF SIZE 2 & 3 FOR NOW', PG, NPERMSIZE(PG)
          STOP
       ENDIF

       START = START + NPERMSIZE(PG)
       ENDDO

       DISTANCE = DOT_PRODUCT(RF-RS,RF-RS)
       IF (PTEST) PRINT*, "MINPERM_CHIRAL> FINAL DISTANCE", DISTANCE


      END SUBROUTINE MINPERM_CHIRAL
