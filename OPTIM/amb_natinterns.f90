! MSB50 --> EFK_ CHARMM/NATINTERNS.SRC


      SUBROUTINE AMBGETNATINTERN
! GENERATE THE LISTS SPECIFYING THE NATURAL INTERNAL COORDINATES FOR THE SYSTEM
      USE COMMONS, ONLY : NATOMS
      USE MODAMBER9
      USE INTCOMMONS
      
      IMPLICIT NONE

      INTEGER A, B, I, R
      CHARACTER*4 PHERING(6), PRORING(5), TRPRING5(5),  &
     &     TRPRING6(6), TRPFRNG(6), HISRING(5), TMPTYPE, &
     &     CRING(6), ARING5(5), ARING6(6), AFRNG(6), RRING(5)
      INTEGER RINGID(NATOMS)   !WHICH RING IS THIS ATOMS PART OF
      INTEGER ATM1, ATM2, RNUM
      INTEGER PATATM, PATRES
      CHARACTER PTYP
      INTEGER NTERM
      CHARACTER(LEN=4) :: RES_RNUM
      LOGICAL ALLRING, RINGATM(NATOMS)
      LOGICAL BBATM1, BBATM2, BBONLY
      INTEGER GLYCA(NATOMS), NINDEX !MSB50
      LOGICAL NOT_NTERM, RINGNEIGHB

      ! TESTING ONLY
      INTEGER TEST

      GLYCA(:)=0
      !THIS ARRAY IS TO SHOW WHICH CA'S ARE GLYCINE I.E. NOT BACKBONE

      !MSB50 - TRY TO FAKE CHARMM
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

!      !MSB50
!      IF(.NOT.ALLOCATED(IBIB)) THEN
!      ALLOCATE(IBIB(NBONA+NBONH))
!      ALLOCATE(JBJB(NBONA+NBONH))
!      IBIB(1) = 1; JBJB(1) = 2
!      IBIB(2) = 2; JBJB(2) = 4
!      IBIB(3) = 2; JBJB(3) = 3
!      IBIB(4) = 4; JBJB(4) = 6
!      IBIB(5) = 4; JBJB(5) = 5
!      ENDIF

      CALL AMB_GETBACKBONE
      IF (.NOT.ALLOCATED(MBBONDNUM)) ALLOCATE(MBBONDNUM(NATOMS))
      IF (.NOT.ALLOCATED(MBADJACENT)) ALLOCATE(MBADJACENT(NATOMS,6)) 
! MSB50 DEALLOCATE IN FINDPERMDIH
      MBBONDNUM(:)=0
      NCNT2 = 0
      NGLYDIH=0     
 
      ! RING DEFINITIONS FOR KNOWN RING RESIDUES
      RRING = (/"C4' ", "O4' ", "C1' ", "C2' ", "C3' "/)
      ARING5 = (/"N9  ", "C8  ", "N7  ", "C5  ", "C4  "/)
      ARING6 = (/"C5  ", "C6  ", "N1  ", "C2  ", "N3  ", "C4  "/)
      AFRNG = (/"C4  ", "C5  ", "N7  ", "N9  ", "N3  ", "C6  "/)
      CRING = (/'N1  ', 'C6  ', 'C5  ', 'C4  ', 'N3  ', 'C2  '/)
      !GRING5 = (/'N9 ', 'C8  ', 'N7  ', 'C5  ', 'C4  '/)
      !GRING6 = (/'C4  ', 'C5  ', 'C6  ', 'N1  ', 'C2  ', 'N3  '/)
      !GFRNG = (/'C4  ', 'C5  ', 'N9  ', 'N7  ', 'N3  ', 'C6  '/)
      !URING = (/'N1  ', 'C6  ', 'C5  ', 'C4  ', 'N3  ', 'C2  '/)
      PHERING = (/'CG  ', 'CD1 ', 'CE1 ', 'CZ  ', 'CE2 ', 'CD2 '/)
      PRORING = (/'N   ', 'CA  ', 'CB  ', 'CG  ', 'CD  '/)
      TRPRING5 = (/'CG  ', 'CD2 ', 'CE2 ', 'NE1 ', 'CD1 '/)
      TRPRING6 = (/'CD2 ', 'CE3 ', 'CZ3 ', 'CH2 ', 'CZ2 ', 'CE2 '/)
      TRPFRNG = (/'CE2 ', 'CD2 ', 'NE1 ', 'CG  ', 'CZ2 ', 'CE3 '/)
      HISRING = (/'CG  ', 'ND1 ', 'CE1 ', 'NE2 ', 'CD2 '/) 

      IF (NUBONDS.GT.0) THEN
         ! ADD ON USER-SPECIFIED BONDS
         DO B = 1,NUBONDS
            IBIB(NBOND+B) = UBONDS(B,1)
            JBJB(NBOND+B) = UBONDS(B,2)
         ENDDO
         NBOND = NBOND + NUBONDS
      ENDIF

      ! LIST NUMBER OF BONDS AND ADJACENT ATOMS FOR EACH ATOM
      NBDS = 0
      DO B = 1,NBOND
         ATM1 = IBIB(B); ATM2 = JBJB(B)

         IF (.NOT.USECART(ATM1).OR..NOT.USECART(ATM2)) NBDS = NBDS + 1

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

      NINTC = NBDS + 3*NCRT

      NRNG = 0; NFRG = 0; NDIH = 0
      RINGATM(:) = .FALSE.
      RINGS(:,:) = -1
      RINGID(:) =-1
      DO RNUM = 1,CARTRESSTART - 1
         ! PULL OUT RINGS FROM THE RING RESIDUES
         ! LABEL RING ATOMS IN RINGATM LIST
         RES_RNUM = IH(M02+RNUM-1)
         IF (RES_RNUM.EQ.'PHE'.OR.RES_RNUM.EQ.'TYR' &
     &      .OR.RES_RNUM.EQ.'NPHE'.OR.RES_RNUM.EQ.'NTYR' &
     &      .OR.RES_RNUM.EQ.'CPHE'.OR.RES_RNUM.EQ.'CTYR') THEN
            NINTC = NINTC + 6
            NDIH = NDIH + 6 ! KEEP TRACK OF INDIVIDUAL DIHEDRALS
            NRNG = NRNG+1
            DO A = 1,6               
               CALL AMB_PATOM(RINGS(NRNG,A), RNUM, PHERING(A))
               IF (RINGS(NRNG,A).LT.0) THEN
                  PRINT*, 'ERROR! COULD NOT FIND ATOM: ', PHERING(A), RNUM
                  STOP
               ENDIF
               RINGATM(RINGS(NRNG,A)) = .TRUE.
               RINGID(RINGS(NRNG,A)) = NRNG
            ENDDO
         ELSE IF (RES_RNUM.EQ.'HIE'.OR.RES_RNUM.EQ.'HID'.OR.  &
     &          RES_RNUM.EQ.'HSP'    &
     &          .OR.RES_RNUM.EQ.'NHIE'.OR.RES_RNUM.EQ.'NHID'.OR. &
     &          RES_RNUM.EQ.'NHSP'.OR.RES_RNUM.EQ.'CHIE'.OR.RES_RNUM &
     &          .EQ.'CHID'.OR. RES_RNUM.EQ.'CHSP') THEN
            NINTC = NINTC + 4
            NRNG = NRNG + 1
            NDIH = NDIH + 5
            DO A = 1,5
               CALL AMB_PATOM(RINGS(NRNG,A), RNUM, HISRING(A))
               IF (RINGS(NRNG,A).LT.0) THEN
                  PRINT*, 'ERROR! COULD NOT FIND ATOM: ', HISRING(A), RNUM
                  STOP
               ENDIF
               RINGATM(RINGS(NRNG,A)) = .TRUE.
               RINGID(RINGS(NRNG,A)) = NRNG
            ENDDO
         ELSE IF ((RES_RNUM.EQ.'PRO'.OR.RES_RNUM.EQ.'NPRO' &
     &        .OR. RES_RNUM.EQ.'CPRO').AND..NOT.BBCART) THEN
            NINTC = NINTC + 4
            NRNG = NRNG + 1 
            NDIH = NDIH + 5
            DO A = 1,5
               CALL AMB_PATOM(RINGS(NRNG,A), RNUM, PRORING(A))
               IF (RINGS(NRNG,A).LT.0) THEN
                  PRINT*, 'ERROR! COULD NOT FIND ATOM: ', PRORING(A), RNUM
                  STOP
               ENDIF
               RINGATM(RINGS(NRNG,A)) = .TRUE.
               RINGID(RINGS(NRNG,A)) = NRNG
            ENDDO
         ELSE IF (RES_RNUM.EQ.'TRP'.OR.RES_RNUM.EQ.'NTRP'&
     &      .OR.RES_RNUM.EQ.'CTRP') THEN
            NINTC = NINTC + 11
            NRNG = NRNG + 1
            NDIH = NDIH + 13
            DO A = 1,5
               CALL AMB_PATOM(RINGS(NRNG,A), RNUM, TRPRING5(A))
               IF (RINGS(NRNG,A).LT.0) THEN
                  PRINT*, 'ERROR! COULD NOT FIND ATOM: ', TRPRING5(A), RNUM
                  STOP
               ENDIF
               RINGATM(RINGS(NRNG,A)) = .TRUE.
               RINGID(RINGS(NRNG,A)) = NRNG
            ENDDO
            NRNG = NRNG + 1 
            DO A = 1,6
               CALL AMB_PATOM(RINGS(NRNG,A), RNUM, TRPRING6(A))
               IF (RINGS(NRNG,A).LT.0) THEN
                  PRINT*, 'ERROR! COULD NOT FIND ATOM: ', TRPRING5(A), RNUM
                  STOP
               ENDIF
               RINGATM(RINGS(NRNG,A)) = .TRUE.
               RINGID(RINGS(NRNG,A)) = NRNG -1 !HAVE SAME RINGID!!!
            ENDDO
            NFRG = NFRG + 1
            DO A = 1,6
               CALL AMB_PATOM(FRINGS(NFRG,A), RNUM, TRPFRNG(A))
               IF (FRINGS(NFRG,A).LT.0) THEN
                  PRINT*, 'ERROR! COULD NOT FIND ATOM: ', TRPRING5(A), RNUM
                  STOP
               ENDIF
            ENDDO
         !MSB50 GLYCINE
         ELSE IF (RES_RNUM.EQ.'GLY'.OR.RES_RNUM.EQ.'NGLY'&
     &       .OR.RES_RNUM.EQ.'CGLY') THEN
             DO A = IX(I02+RNUM-1), IX(I02+RNUM)-1
                IF (IH(M04+A-1).EQ.'CA  ') THEN
                  GLYCA(NGLYDIH+1)=A
                  NGLYDIH=NGLYDIH+1
                ENDIF
             ENDDO
!MSB50 RNA
         ELSE IF (RES_RNUM.EQ.'RA  '.OR.RES_RNUM.EQ.'RA3 '&
             .OR.RES_RNUM.EQ.'RA5 '.OR.RES_RNUM.EQ.'RAN ' &
             .OR. RES_RNUM.EQ.'RG  '.OR.RES_RNUM.EQ.'RG3 ' &
             .OR.RES_RNUM.EQ.'RG5 '.OR.RES_RNUM.EQ.'RGN ') THEN
             NINTC = NINTC +11
             NRNG = NRNG + 1
             NDIH = NDIH + 13
             DO A = 1,5
                CALL AMB_PATOM(RINGS(NRNG,A), RNUM, ARING5(A))
                IF (RINGS(NRNG, A).LT.0)  THEN
                   PRINT*, 'ERROR! COULD NOT FIND ATOM: ', ARING5(A), RNUM
                   STOP
                ENDIF
                RINGATM(RINGS(NRNG,A)) =.TRUE.
                RINGID(RINGS(NRNG,A)) = NRNG
             ENDDO
             NRNG = NRNG + 1
             DO A =1,6
                CALL AMB_PATOM(RINGS(NRNG, A), RNUM, ARING6(A))
                IF (RINGS(NRNG, A).LT.0) THEN
                   PRINT*, 'ERROR! COULD NOT FIND ATOM: ', ARING6(A), RNUM
                   STOP
                ENDIF
                RINGATM(RINGS(NRNG,A)) =.TRUE.
                RINGID(RINGS(NRNG,A)) = NRNG - 1 ! HAVE SAME RINGID AS SAME RING
             ENDDO
             NFRG = NFRG + 1
             DO A = 1,6
               CALL AMB_PATOM(FRINGS(NFRG,A), RNUM,  AFRNG(A))
               IF (FRINGS(NFRG,A).LT.0) THEN
                  PRINT*, 'ERROR! COULD NOT FIND ATOM: ', AFRNG(A),RNUM
                  STOP
               ENDIF
             ENDDO
!RIBOSE
             NINTC = NINTC + 4
             NRNG = NRNG + 1
             NDIH = NDIH + 5
             DO A = 1,5
               CALL AMB_PATOM(RINGS(NRNG,A), RNUM, RRING(A))
               IF (RINGS(NRNG,A).LT.0) THEN
                  PRINT*, 'ERROR! COULD NOT FIND ATOM: ', RRING(A), RNUM
                  STOP
               ENDIF
               RINGATM(RINGS(NRNG,A)) = .TRUE.
               RINGID(RINGS(NRNG,A)) = NRNG
            ENDDO
         ELSE IF (RES_RNUM.EQ.'RC  '.OR.RES_RNUM.EQ.'RC3 '&
             .OR.RES_RNUM.EQ.'RC5 '.OR.RES_RNUM.EQ.'RCN ' &
             .OR. RES_RNUM.EQ.'RU  '.OR.RES_RNUM.EQ.'RU3 ' &
             .OR.RES_RNUM.EQ.'RU5 '.OR.RES_RNUM.EQ.'RUN ') THEN
            NINTC = NINTC + 6
            NDIH = NDIH + 6 ! KEEP TRACK OF INDIVIDUAL DIHEDRALS
            NRNG = NRNG+1
            DO A = 1,6
               CALL AMB_PATOM(RINGS(NRNG,A), RNUM, CRING(A))
               IF (RINGS(NRNG,A).LT.0) THEN
                  PRINT*, 'ERROR! COULD NOT FIND ATOM: ', CRING(A), RNUM
                  STOP
               ENDIF
               RINGATM(RINGS(NRNG,A)) = .TRUE.
               RINGID(RINGS(NRNG,A)) = NRNG
            ENDDO
!RIBOSE
             NINTC = NINTC + 4
             NRNG = NRNG + 1
             NDIH = NDIH + 5
             DO A = 1,5
               CALL AMB_PATOM(RINGS(NRNG,A), RNUM, RRING(A))
               IF (RINGS(NRNG,A).LT.0) THEN
                  PRINT*, 'ERROR! COULD NOT FIND ATOM: ', RRING(A), RNUM
                  STOP
               ENDIF
               RINGATM(RINGS(NRNG,A)) = .TRUE.
               RINGID(RINGS(NRNG,A)) = NRNG
             ENDDO

         ENDIF
      ENDDO
 
      !STORE ANGLES FOR GLYCINE
      IF (.NOT.ALLOCATED(GLYDIH)) ALLOCATE(GLYDIH(NGLYDIH,8))
      DO I=1,NGLYDIH
        NINDEX = 0
        NOT_NTERM = .FALSE.
        GLYDIH(I,1) = GLYCA(I)
        A = GLYCA(I)  
        DO B=1,MBBONDNUM(GLYCA(I)) 
           IF (MBBONDNUM(MBADJACENT(A,B)).EQ.1) THEN
              GLYDIH(I,8-NINDEX) = MBADJACENT(A,B)
              NINDEX = NINDEX +1
           ELSEIF (IH(M04+MBADJACENT(A,B)-1).EQ.'C   ') THEN
              GLYDIH(I,2) = MBADJACENT(A,B)
              GLYDIH(I,4) = MBADJACENT(MBADJACENT(A,B),1)
              GLYDIH(I,5) = MBADJACENT(MBADJACENT(A,B),2)
           ELSE !N
              GLYDIH(I,3) = MBADJACENT(A,B)
              DO R = 1,MBBONDNUM(MBADJACENT(A,B))
                IF (MBBONDNUM(MBADJACENT(MBADJACENT(A,B),R)).GT.1&
     &           .AND.MBADJACENT(MBADJACENT(A,B),R).NE.A) THEN 
                   GLYDIH(I,6) = MBADJACENT(MBADJACENT(A,B),R)
                   NOT_NTERM = .TRUE.
                ENDIF
              ENDDO
              IF (.NOT.NOT_NTERM) THEN
                 IF (MBADJACENT(MBADJACENT(A,B),1).NE.A) THEN
                    GLYDIH(I,6) = MBADJACENT(MBADJACENT(A,B),1)
                 ELSE
                    GLYDIH(I,6) = MBADJACENT(MBADJACENT(A,B),2)
                 ENDIF
              ENDIF
            ENDIF
         ENDDO
      ENDDO
      

      ! PUT IN THE EXTRA (USER-SPECIFIED) RINGS
      DO R = 1,NURINGS
         NRNG = NRNG + 1
         NDIH = NDIH + URINGS(R,0)
         IF (URINGS(R,0).EQ.5) THEN
            NINTC = NINTC + 4
         ELSE IF (URINGS(R,0).EQ.6) THEN
            NINTC = NINTC + 6
         ELSE
            PRINT*, 'GETNATINTERNS>> ERROR!&
              & USER-SPECIFIED RINGS MUST HAVE 5 OR 6 ATOMS'
            STOP
         ENDIF
            
         DO I = 1, URINGS(R,0)
            RINGS(NRNG,I) = URINGS(R,I)
            RINGATM(URINGS(R,I)) = .TRUE.
         ENDDO
      ENDDO

      ! PUT IN THE NON-TERMINAL CENTERS
      RNUM = 0
      NCNT = 0; NIMP = 0
      DO A = 1,CARTATMSTART - 1
         IF (RNUM.LT.NRES) THEN
            IF (A.EQ.IX(I02+RNUM)) THEN
                RNUM = RNUM + 1 !NEW RESIDUE STARTED
             ENDIF
         ENDIF
         RES_RNUM = IH(M02+RNUM-1)        

         IF (MBBONDNUM(A).EQ.0) THEN
            PRINT*, 'GETNATINT>> ERROR! DISCONNECTED ATOM.', A
            STOP
         ELSE IF (MBBONDNUM(A).EQ.1) THEN
            CYCLE ! IGNORE TERMINAL ATOMS
         ENDIF

         IF (RINGATM(A)) THEN
            ! IGNORE CENTERS WHERE ALL SUBSTITUENTS ARE RING ATMS
            ALLRING = .TRUE.
            DO I = 1,MBBONDNUM(A)
               IF (.NOT.RINGATM(MBADJACENT(A,I)) .OR. &
     &            (RINGID(A).NE.RINGID(MBADJACENT(A,I)))) THEN
                  ALLRING = .FALSE.              
                  EXIT
               ENDIF
            ENDDO
            IF (ALLRING) THEN
               CYCLE 
            ENDIF
         ENDIF

         NCNT = NCNT + 1
         CENTERS(NCNT,1) = A

         NTERM = 0 ! NUMBER OF TERMINAL OR QUASI-TERMINAL NEIGHBORS        
         RINGNEIGHB = .FALSE.
         DO I=1, MBBONDNUM(A)
            IF (RINGATM(A).AND.(RINGID(MBADJACENT(A,I)).LT.RINGID(A))) THEN
               RINGNEIGHB=.TRUE. 
               CYCLE
            ENDIF
         ENDDO
         !PREVIOUS LOOP ONLY DETERMINED WHETHER TWO RINGS ARE LINKED

         DO I = 1,MBBONDNUM(A)
            ! TREAT ALL RING SUBSTITUENTS AS TERMINAL
            ! CB IS QUASITERMINAL FOR CA CENTER
            ! THE NH'S OF ARG ARE QUASI-TERMINAL
            !MSB50 MAKE O5 QUASITERMINAL FOR P TO MAKE P METHINE SP3
            ! LINKED RINGS: LINKATOM IN FOLLOWING RING QUASI-TERMINAL
            ! ITS NEIGHBOURS ALSO QUASITERMINAL TO FIX FOLLOWING RING POSITION
     !&           .OR.(RINGATM(A).AND.RINGNEIGHB.AND.                    &
     !&              (RINGID(MBADJACENT(A,I)).EQ.RINGID(A))).OR.         &
            IF (MBBONDNUM(MBADJACENT(A,I)).EQ.1.OR.                     & 
     &           (RINGATM(A).AND..NOT.RINGATM(MBADJACENT(A,I))).OR.     &
     &           (RINGATM(A).AND.(RINGID(MBADJACENT(A,I)).NE.RINGID(A))).OR.&
     &           IH(M04+A-1).EQ.'CA'.AND.IH(M04+(MBADJACENT(A,I))-1).EQ.&
     &           'CB'.OR.                                               &
     &           (RES_RNUM.EQ.'ARG'.AND.                                &
     &           (IH(M04+(MBADJACENT(A,I))-1).EQ.'NH1'.OR.              &
     &           IH(M04+(MBADJACENT(A,I))-1).EQ.'NH2')).OR.             &
     &           (IH(M04+A-1).EQ.'P   '.AND.                            &
     &           IH(M04+(MBADJACENT(A,I))-1).EQ."O5' ")) THEN
               
               NTERM = NTERM + 1
               ! PUT TERMINAL ATMS AT END
               CENTERS(NCNT,MBBONDNUM(A)-NTERM+2) = MBADJACENT(A,I)
               !MSB50 - PROLINE EXCEPTION
               IF ((RES_RNUM.EQ.'PRO'.OR.RES_RNUM.EQ.'NPRO'.OR.         &
     &         RES_RNUM.EQ.'CPRO').AND.IH(M04+MBADJACENT(A,I)-1).EQ.'HA' &
     &             .AND.(MBBONDNUM(A)-NTERM+2).EQ.3) THEN
                    CENTERS(NCNT, 3)= CENTERS(NCNT, 4)
                    CENTERS(NCNT, 4) = MBADJACENT(A,I)
               ENDIF
            ELSE
               CENTERS(NCNT,I-NTERM+1) = MBADJACENT(A,I)
            ENDIF
         ENDDO
          

         ! SET CENTER TYPE
         IF (RINGATM(A)) THEN
            IF (MBBONDNUM(A).EQ.4) THEN
              ! IF (RINGNEIGHB.AND.NTERM == 2 ) THEN
              !    CENTERS(NCNT,0) = 2
              !    NINTC = NINTC + 5 
              ! ELSE IF (RINGNEIGHB.AND.NTERM == 3) THEN
              !    CENTERS(NCNT,0) = 1
              !    NINTC = NINTC + 5 
              ! ELSE !NORMAL - NO LINKED RINGS
                  CENTERS(NCNT,0) = 8
                  NINTC = NINTC + 4
              ! ENDIF
            ELSE IF (MBBONDNUM(A).EQ.3) THEN
               CENTERS(NCNT,0) = 7
               NINTC = NINTC + 2
               NIMP = NIMP + 1
               IMPDIH(NIMP,1) = CENTERS(NCNT,1)
               IMPDIH(NIMP,2) = CENTERS(NCNT,4)
               IMPDIH(NIMP,3:4) = CENTERS(NCNT,2:3)
            ELSE
               PRINT*, 'ERROR! CAN HAVE AT MOST TWO RING SUBSTITUENTS'
               PRINT*, A, MBBONDNUM(A)
               STOP
            ENDIF
         ELSE IF (MBBONDNUM(A).EQ.4) THEN
            NINTC = NINTC + 5
            SELECT CASE(NTERM)
            CASE (3)
               CENTERS(NCNT,0) = 1
            CASE (2)
               CENTERS(NCNT,0) = 2
               IF (INTMINPERMT) THEN !MSB50
                   IF (.NOT.ALLOCATED(CENTER2)) THEN
                       ALLOCATE(CENTER2(NATOMS))
                       CENTER2(:) = 0
                   ENDIF
                   IF (IH(M04+A-1).NE.'CA  ') THEN!NO CAS, ONLY IN AM ACID SIDECHAINS
                      NCNT2 = NCNT2+1
                      CENTER2(NCNT2) = NCNT 
                   ENDIF
               ENDIF
            CASE (1)
               CENTERS(NCNT,0) = 3
            CASE DEFAULT
               PRINT*, 'ERROR! UNKNOWN CENTER TYPE. ATM,BONDNUM,NTERM: '&
     &             ,A, MBBONDNUM(A), NTERM
               STOP
            END SELECT
         ELSE IF (MBBONDNUM(A).EQ.3) THEN           
            NIMP = NIMP + 1
            NINTC = NINTC + 3
            SELECT CASE(NTERM)
            CASE (2)
               CENTERS(NCNT,0) = 4
               IMPDIH(NIMP,1:4) = CENTERS(NCNT,1:4)
            CASE (1)
               CENTERS(NCNT,0) = 5
               IMPDIH(NIMP,1) = CENTERS(NCNT,1)
               IMPDIH(NIMP,2) = CENTERS(NCNT,4)
               IMPDIH(NIMP,3:4) = CENTERS(NCNT,2:3)
            CASE DEFAULT
               PRINT*, 'ERROR! UNKNOWN CENTER TYPE. ATM,BONDNUM,NTERM: ', &
     &              A, MBBONDNUM(A), NTERM
               STOP
            END SELECT
         ELSE IF (MBBONDNUM(A).EQ.2) THEN
            NINTC = NINTC + 1
            CENTERS(NCNT,0) = 6
         ENDIF

         ! COUNTS AS BACKBONE CENTER ONLY IF FIRST 3 ATOMS ARE IN BACKBONE
         IF (BBCART) THEN
            BBONLY = .TRUE.
            DO I = 1,3
               IF (.NOT.AM_BACKBONE(CENTERS(NCNT,I))) BBONLY = .FALSE.
            ENDDO
           
            IF (BBONLY) THEN
               SELECT CASE (CENTERS(NCNT,0))
               CASE(5)          !TREAT LIKE RING SUBSTITUENT
                  CENTERS(NCNT,0) = 7
                  NINTC = NINTC - 1
               CASE(2)          ! TREAT LIKE RING SUBSTITUENT
                  CENTERS(NCNT,0) = 8
                  NINTC = NINTC - 1
               CASE(6)          !SKIP COORD
                  NINTC = NINTC -1
                  NCNT = NCNT - 1
               CASE(7)          !SKIP COORD
                  NINTC = NINTC - 2
                  NCNT = NCNT - 1
               CASE DEFAULT
                  PRINT*, 'ERROR! IF BBCART SET,'
                  PRINT*, 'EXPECT ALL BACKBONE CENTERS TO BE'
                  PRINT*, 'OF TYPE 2,5,6,OR 7', RNUM, RES_RNUM
                  PRINT*, CENTERS(NCNT,:)
                  PRINT*, "MSB50 - AND 8 - FIX THIS IF IT'S 8"
                  STOP
               END SELECT
            ENDIF
         ENDIF
      ENDDO

      ! GET LINEAR DIHEDRALS
      NLDH = 0
      DO B = 1,NBOND
         ATM1 = IBIB(B); ATM2 = JBJB(B)
 
         ! IGNORE BONDS INVOLVING A TERMINAL ATOM, BONDS IN RINGS
         IF ((MBBONDNUM(ATM1).EQ.1.OR.MBBONDNUM(ATM2).EQ.1).OR.&
     &        (RINGATM(ATM1).AND.RINGATM(ATM2).AND.&
     &         RINGID(ATM1).EQ.RINGID(ATM2))) CYCLE        
             !LINDH REQUIRED IF LINKED RINGS 

         NLDH = NLDH + 1
         LINDIH(NLDH,1) = ATM1; LINDIH(NLDH,2) = ATM2


         BBATM1 = IH(M04+ATM1-1).EQ.'N'.OR.IH(M04+ATM1-1).EQ.'CA'.OR.   & 
     &       IH(M04+ATM1-1).EQ.'C'.OR.IH(M04+ATM1-1).EQ.'CH3 '.OR.      &
     &       IH(M04+ATM1-1).EQ."C5' ".OR.IH(M04+ATM1-1).EQ."C3' ".OR.   &
     &       IH(M04+ATM1-1).EQ."O3' ".OR.IH(M04+ATM1-1).EQ."P"
         BBATM2 = IH(M04+ATM2-1).EQ.'N'.OR.IH(M04+ATM2-1).EQ.'CA'.OR.   &
     &        IH(M04+ATM2-1).EQ.'C'.OR.IH(M04+ATM2-1).EQ.'CH3 '.OR.     &
     &       IH(M04+ATM2-1).EQ."C5' ".OR.IH(M04+ATM2-1).EQ."C3' ".OR.   &
     &       IH(M04+ATM2-1).EQ."O3' ".OR.IH(M04+ATM2-1).EQ."P"
         BBONLY = .FALSE.


         ! WHEN APPLICABLE, USE BACKBONE BONDS ONLY
         IF (BBATM1.AND.BBATM2) THEN
            BBONLY = .FALSE.
            DO I = 1,MBBONDNUM(ATM1)
               TMPTYPE = IH(M04+(MBADJACENT(ATM1,I))-1)
               IF (TMPTYPE.NE.IH(M04+ATM2-1).AND.                       &
     &              (TMPTYPE.EQ.'N'.OR.TMPTYPE.EQ.'CA'.OR.              &
     &              TMPTYPE.EQ.'C'.OR.TMPTYPE.EQ.'CH3 '.OR.             &
     &              TMPTYPE.EQ."C5' ".OR.TMPTYPE.EQ."C3' ".OR.           &
     &              TMPTYPE.EQ."O3' ".OR.TMPTYPE.EQ."P   ".OR.          &
     &              TMPTYPE.EQ."O5' ")) THEN
                  LINDIH(NLDH,-1) = 1
                  LINDIH(NLDH,3) = MBADJACENT(ATM1,I)
                  BBONLY = .TRUE.
                  EXIT
               ENDIF
            ENDDO
         ENDIF
         IF (.NOT.BBONLY) THEN
            LINDIH(NLDH,-1)= MBBONDNUM(ATM1) - 1
            A = 2
            DO I = 1,MBBONDNUM(ATM1)
               IF (MBADJACENT(ATM1,I).NE.ATM2) THEN
                  A = A + 1
                  LINDIH(NLDH,A) = MBADJACENT(ATM1,I)
               ENDIF
            ENDDO
         ENDIF

         IF (BBATM1.AND.BBATM2) THEN
            BBONLY = .FALSE.
            DO I = 1,MBBONDNUM(ATM2)
               TMPTYPE = IH(M04+MBADJACENT(ATM2,I)-1)
               IF (TMPTYPE.NE.IH(M04+ATM1-1).AND.                       &
     &              (TMPTYPE.EQ.'N'.OR.TMPTYPE.EQ.'CA'.OR.              &
     &                TMPTYPE.EQ.'C'.OR.TMPTYPE.EQ.'CH3 '.OR.           &
     &              TMPTYPE.EQ."C5' ".OR.TMPTYPE.EQ."C3' ".OR.          &
     &              TMPTYPE.EQ."O3' ".OR.TMPTYPE.EQ."P   ".OR.          &
     &              TMPTYPE.EQ."O5' ")) THEN 
                  LINDIH(NLDH,0) = 1
                  LINDIH(NLDH,LINDIH(NLDH,-1)+3) = MBADJACENT(ATM2,I)
                  BBONLY = .TRUE.
                  EXIT
               ENDIF
            ENDDO
         ENDIF
         IF (.NOT.BBONLY) THEN
            LINDIH(NLDH,0)= MBBONDNUM(ATM2) - 1
            A = LINDIH(NLDH,-1)+2
            DO I = 1,MBBONDNUM(ATM2)
               IF (MBADJACENT(ATM2,I).NE.ATM1) THEN
                  A = A + 1
                  LINDIH(NLDH,A) = MBADJACENT(ATM2,I)
               ENDIF
            ENDDO
         ENDIF

         ! IGNORE IF USING CARTESIANS FOR ALL ATOMS INVOLVED
         BBONLY = .TRUE.
         DO I = 1,LINDIH(NLDH,-1)+LINDIH(NLDH,0)+2
            IF (.NOT.USECART(LINDIH(NLDH,I))) THEN
               BBONLY = .FALSE.
            ENDIF
         ENDDO
      
         IF (BBONLY) THEN
            NLDH = NLDH - 1
         ELSE
            NDIH = NDIH + LINDIH(NLDH,-1)*LINDIH(NLDH,0)
            NINTC = NINTC + 1
         ENDIF
                  
      ENDDO

      IF (PRINTCOORDS) THEN
         PRINT*, 'TOTAL COORDS:', NINTC
         PRINT*, 'CARTESIAN ATOMS:', NCRT

         DO I = 1, NCRT
            PRINT*, CARTATMS(I)
         ENDDO

         PRINT*, 'BONDS:', NBDS
         DO I=1,NBOND
            IF (.NOT.(USECART(IBIB(I)).AND.USECART(JBJB(I))))            &
     &            PRINT*, IBIB(I), JBJB(I)
         ENDDO
         PRINT*, 'CENTERS:', NCNT
         DO I=1,NCNT
            PRINT '(I6,A,6I6)', I, '|', CENTERS(I,0:5)
         ENDDO
         PRINT*, 'RINGS:', NRNG
         DO I=1,NRNG
            PRINT*, RINGS(I,1:6)
            PRINT*, '---'
         ENDDO
         PRINT*, 'FUSED RINGS:', NFRG
         DO I=1,NFRG
            PRINT*, FRINGS(I,1:6)
            PRINT*, '---'
         ENDDO
         PRINT*, 'LINEAR DIHEDRALS', NLDH, NDIH
         DO I=1,NLDH
            PRINT*, I, '|', LINDIH(I,-1:8)
            PRINT*, '---'
         ENDDO
         PRINT*, 'IMPROPER DIHEDRALS', NIMP
         DO I=1,NIMP
            PRINT*, IMPDIH(I,:)
         ENDDO
      ENDIF

      IF ((BBCART.AND.NINTC.NE.3*NATOMS).OR.                            &
     &          (.NOT.BBCART.AND.NINTC.NE.3*NATOMS-6)) THEN
         PRINT*, 'ERROR!'
         PRINT*, 'NUMBER OF NATURAL INTERNALS SHOULD BE 3*NATOMS-6,     &
     &        OR 3*NATOMS IF BBCART', BBCART, NINTC, NATOMS, 3*NATOMS-6
      ENDIF

      RETURN
      END SUBROUTINE


!MSB50 - FOR AMBER
      SUBROUTINE AMB_GETBACKBONE
       ! WHICH OF THE ATOMS ARE PART OF THE BACKBONE?
      USE MODAMBER9 
      USE COMMONS, ONLY: NATOMS
      USE INTCOMMONS, ONLY: INTMINPERMT
 
         IMPLICIT NONE
 
         INTEGER :: A, R
       
         IF(.NOT.ALLOCATED(AM_BACKBONE)) ALLOCATE(AM_BACKBONE(NATOMS))
         AM_BACKBONE(1:NATOMS) = .FALSE.
       DO A = 1,NATOMS
          IF (IH(M04+A-1).EQ.'C   '.OR.IH(M04+A-1).EQ.'CA  '.OR.    &
     &       IH(M04+A-1).EQ.'N   '.OR.IH(M04+A-1).EQ.'CH3 '.OR. &
     &       IH(M04+A-1).EQ."C5' ".OR.IH(M04+A-1).EQ."C3' ".OR. &
     &       IH(M04+A-1).EQ."O3' ".OR.IH(M04+A-1).EQ."P") AM_BACKBONE(A) = .TRUE.
             !LAST CASE - TERMINAL - NEED TO DOUBLE CHECK THIS
       ENDDO

       IF ((IH(M04+A-1).EQ."C5' ".OR.IH(M04+A-1).EQ."C3' ".OR. &
     &       IH(M04+A-1).EQ."O3' ".OR.IH(M04+A-1).EQ."P").AND. &
     &       INTMINPERMT) THEN
           PRINT*, "*************ERROR*******************"
           PRINT*,"INTMINPERM DOES NOT WORK WITH RNA BACKBONE SETTINGS!"
           PRINT*,"USE NOPERMPROCHIRAL IN ADDITION!!"
        ENDIF 

       ! LABEL ALL PROLINE ATOMS (EXCEPT O) AS BACKBONE
       DO R = 1,NRES
          IF ((IH(M02+R-1)).EQ.'PRO'.OR. (IH(M02+R-1)).EQ.'NPRO'.OR.     &
     &        (IH(M02+R-1)).EQ.'CPRO') THEN
             PRINT*, "YAY PROLINE"
             DO A = IX(I02+R-1),IX(I02+R)-1
                IF(IH(M04+A-1)(1:1).NE.'O') AM_BACKBONE(A) = .TRUE.
             ENDDO
          ENDIF
       ENDDO
       
       RETURN
 
      END SUBROUTINE AMB_GETBACKBONE

! ***********************************************************************
!YAY FOR AMBER - AGAIN - SEE EFK
      SUBROUTINE AMB_NATINTSETUP
        ! SETUP UP VARIOUS ARRAYS NECESSARY FOR NATURAL INTERNALS
        ! THIS IS THE ROUTINE THAT DEFINES THE NATURAL INTERNAL COORDINATES IN TERMS OF ANGLE AND TORSION COEFFICIENTS
         USE MODAMBER9
         USE INTCOMMONS
         USE COMMONS, ONLY:NATOMS   
  
         IMPLICIT NONE
         INTEGER, PARAMETER :: LARGENEG = -99999999
         DOUBLE PRECISION :: S6R, S2R, S26R, S18R, S12R, A, B, SABR, SABR2
         
         NBDS = 0; NCNT = 0; NRNG = 0; NFRG = 0; NLDH = 0; NIMP = 0
         NCRT = 0; NDIH = 0
         
         !MSB50 RINGS(3*NRES,6) AS CAN HAVE RIBOSE AND 2 JOINED RINGS PER RESIDUE
         ALLOCATE(CENTERS(NATOMS,0:5), RINGS(3*NRES,6), FRINGS(NRES,6))
         ALLOCATE(LINDIH(NATOMS,-1:8), IMPDIH(NATOMS,4), CARTATMS(NATOMS))
         
         IF (CARTRESSTART.EQ.0) THEN
             PRINT*, "AMB_NATINTSETUP> START FROM RESIDUE 1, NOT 0!"
             CARTRESSTART=1
         ENDIF

         IF (CARTRESSTART.GT.NRES+1.OR.CARTRESSTART.LT.0) THEN
            CARTRESSTART = NRES + 1
            CARTATMSTART = NATOMS + 1
         ELSE
            CARTATMSTART = IX(I02+CARTRESSTART-1)
         ENDIF
         
         !     SET UP COEFF TABLE
         S6R = 1/SQRT(6.0)
         S2R = 1/SQRT(2.0)
         S26R = 1/SQRT(26.0)
         S18R = 1/SQRT(18.0)
         S12R = 1/SQRT(12.0)

         COEFF(1,1,1:6) = (/-1,-1,-1,1,1,1/)*S6R
         COEFF(1,2,1:6) = (/0,0,0,2,-1,-1/)*S6R
         COEFF(1,3,1:6) = (/0,0,0,0,1,-1/)*S2R
         COEFF(1,4,1:6) = (/2,-1,-1,0,0,0/)*S6R
         COEFF(1,5,1:6) = (/0,1,-1,0,0,0/)*S2R     
       
         COEFF(2,1,1:6) = (/0,1,-1,1,-1,0/)*0.5
         COEFF(2,2,1:6) = (/0,-1,-1,1,1,0/)*0.5
         COEFF(2,3,1:6) = (/0,-1,1,1,-1,0/)*0.5
         COEFF(2,4,1:6) = (/1,0,0,0,0,5/)*S26R
         COEFF(2,5,1:6) = (/5,0,0,0,0,1/)*S26R
       
         COEFF(3,1,1:6) = (/0,0,2,0,-1,-1/)*S6R
         COEFF(3,2,1:6) = (/0,0,0,0,1,-1/)*S2R
         COEFF(3,3,1:6) = (/1,4,0,1,0,0/)*S18R
         COEFF(3,4,1:6) = (/1,1,0,4,0,0/)*S18R
         COEFF(3,5,1:6) = (/4,1,0,1,0,0/)*S18R
       
         COEFF(4,1,1:3)=(/-1,-1,2/)*S6R
         COEFF(4,2,1:3)=(/1,-1,0/)*S2R
       
         COEFF(5,1,1:3)=(/0,1,-1/)*S2R
         COEFF(5,2,1:3)=(/2,-1,-1/)*S6R
       
         COEFF(6,1,1:1) = (/1/)
       
         COEFF(7,1,1:2) = (/1,-1/)*S2R
       
         COEFF(8,1,1:5) = (/1,-1,1,-1,0/)*0.5
         COEFF(8,2,1:5) = (/-1,-1,1,1,0/)*0.5
         COEFF(8,3,1:5) = (/-1,1,1,-1,0/)*0.5
         COEFF(8,4,1:5) = (/0,0,0,0,1/)
       
         A = COS(144.0*PI/180.0)
         B = COS(72.0*PI/180.0)
         SABR = 1/SQRT(1+2*A*A+2*B*B)
         SABR2 = 1/SQRT(2*(A-B)**2+2*(1-A)**2)
         COEFF(25,1,1:5) = (/1.0D0,A,B,B,A/)*SABR
         COEFF(25,2,1:5) = (/0.0D0,A-B,1.0D0-A,A-1.0D0,B-A/)*SABR2
         COEFF(25,3,1:5) = (/B,A,1.0D0,A,B/)*SABR
         COEFF(25,4,1:5) = (/A-1.0D0,B-A,0.0D0,A-B,1.0D0-A/)*SABR2
       
         COEFF(26,1,1:6) = (/1,-1,1,-1,1,-1/)*S6R
         COEFF(26,2,1:6) = (/2,-1,-1,2,-1,-1/)*S6R
         COEFF(26,3,1:6) = (/0,1,-1,0,1,-1/)*0.5D0
         COEFF(26,4,1:6) = (/1,-1,1,-1,1,-1/)*S6R
         COEFF(26,5,1:6) = (/1,0,-1,1,0,-1/)*0.5D0
         COEFF(26,6,1:6) = (/-1,2,-1,-1,2,-1/)*S12R
       
         !     SET UP ATOMXPOSINC TABLE
         !     LARGE NEGATIVE NUMBER LISTED FOR THOSE ATOMS NOT INVOLVED
         ATOMXPOSINC(1,1,1:5) = (/0,3,6,9,12/)
         ATOMXPOSINC(1,2,1:5) = (/15,LARGENEG,18,21,24/)
         ATOMXPOSINC(1,3,1:5) = (/27,LARGENEG,30,33,36/)
         ATOMXPOSINC(1,4,1:5) = (/39,42,45,48,51/)
         ATOMXPOSINC(1,5,1:5) = (/54,57,LARGENEG,60,63/)

         ATOMXPOSINC(2,1,1:5) = (/0,3,6,9,12/)
         ATOMXPOSINC(2,2,1:5) = (/15,18,21,24,27/)
         ATOMXPOSINC(2,3,1:5) = (/30,33,36,39,42/)
         ATOMXPOSINC(2,4,1:5) = (/45,48,51,54,57/)
         ATOMXPOSINC(2,5,1:5) = (/60,63,66,69,72/)
       
         ATOMXPOSINC(3,1,1:5) = (/0,3,6,9,12/)
         ATOMXPOSINC(3,2,1:5) = (/15,LARGENEG,18,21,24/)
         ATOMXPOSINC(3,3,1:5) = (/27,30,33,36,LARGENEG/)
         ATOMXPOSINC(3,4,1:5) = (/39,42,45,48,LARGENEG/)
         ATOMXPOSINC(3,5,1:5) = (/51,54,57,60,LARGENEG/)
       
         ATOMXPOSINC(4,1,1:4)=(/0,3,6,9/)
         ATOMXPOSINC(4,2,1:4)=(/12,15,18,21/)
       
         ATOMXPOSINC(5,1,1:4)=(/0,3,6,9/)
         ATOMXPOSINC(5,2,1:4)=(/12,15,18,21/)
       
         ATOMXPOSINC(6,1,1:3) = (/0,3,6/)
       
         ATOMXPOSINC(7,1,1:4) = (/0,3,6,9/)      
      
         ATOMXPOSINC(8,1,1:5)=(/0,3,6,9,12/)
         ATOMXPOSINC(8,2,1:5)=(/15,18,21,24,27/)
         ATOMXPOSINC(8,3,1:5)=(/30,33,36,39,42/)
         ATOMXPOSINC(8,4,1:5)=(/45, LARGENEG, LARGENEG,48,51/)
       END SUBROUTINE AMB_NATINTSETUP
       


! ************************************************************************************       
      SUBROUTINE AMB_GETBEENAT(XCART,XINT,VAL,INDX,JNDX,NINTC,NCART,NNZ,GCS,NOCOOR,NODERV,KD)
      USE COMMONS, ONLY : NATOMS
      USE MODAMBER9, ONLY: IBIB,JBJB, NBONA,NBONH,IH,M04
      USE INTCOMMONS, ONLY : CARTATMSTART, PRINTCOORDS, COEFF, CENTERS, RINGS,&
     &     FRINGS, LINDIH, IMPDIH, CARTATMS, ATOMXPOSINC, NBDS, NCNT, NRNG,&
     &     NFRG, NLDH, NIMP, NCRT, USECART, PREVDIH
      
      IMPLICIT NONE
! GET INTERNAL COORDINATES AND DERIVATIVE MATRIX FOR NATURAL INTERNAL COORDS


!-----------------------------------------------------------------------
!
! NNZ IS NUMBER OF NON-ZERO ELEMENTS IN B-MATRIX
! NINTC IS NUMBER OF INTERNAL COORDINATES
! NCART IS NUMBER OF CARTESIANS
!
       INTEGER NNZ,NINTC,NCART
!
       INTEGER INDX(NNZ),JNDX(NNZ)
       INTEGER MMM,ITHH,IPHII,II,JJ,KK,LL,IIKD,JJKD,KKKD,LLKD,KD1,KD
       DOUBLE PRECISION XCART(NCART),XINT(NINTC)
       DOUBLE PRECISION X(NCART),Y(NCART),Z(NCART)
       DOUBLE PRECISION VAL(NNZ),GCS(2*KD+1,NCART)
       DOUBLE PRECISION RX,RY,RZ,S2,S,SR,RXSR,RYSR,RZSR
       DOUBLE PRECISION DXI,DYI,DZI,DXJ,DYJ,DZJ,RI2,RJ2,RIR,RJR,RI,RJ
       DOUBLE PRECISION DXIR,DYIR,DZIR,DXJR,DYJR,DZJR,CST,ISNT
       DOUBLE PRECISION FX,FY,FZ,GX,GY,GZ,HX,HY,HZ,AX,AY,AZ,BX,BY,BZ
       DOUBLE PRECISION CSTTWO,SNTTWO2,SNTTWO2R,CSTTHREE,SNTTHREE2,SNTTHREE2R,DUMMY,DUMMY2
       DOUBLE PRECISION VXI,VYI,VZI,VXJ,VYJ,VZJ,VXK,VYK,VZK,VXL,VYL,VZL
       INTEGER MM,ITH,IPHI,I,J,K,L,IC,I1,J1,SUMPHI,DIFF,ATOMWIDTH
       INTEGER II1,JJ1,KK1,LL1,IIKD1,JJKD1,KKKD1,LLKD1
       INTEGER II2,JJ2,KK2,LL2,IIKD2,JJKD2,KKKD2,LLKD2
       LOGICAL NOCOOR, NODERV
!JMC
       DOUBLE PRECISION MYSCALAR,MYTX,MYTY,MYTZ

      INTEGER CURCOORD, CURBIND, DIHC
      INTEGER XIND(8)
      INTEGER IX, JX, KX, LX
      INTEGER K1, L1, X1, X2, A1, C1, TCOUNT, NANGC, NDIHC
      INTEGER CNT, LDH, IMP, RNG, FRG, CNTATMS, RATMS
      DOUBLE PRECISION CURCOEFF
      INTEGER INFO(-1:8), CNTTYPE, RTYPE
      INTEGER ATOMXSHIFT, ISHIFT, JSHIFT
      DOUBLE PRECISION PHI, THETA, DTDI(0:2), DTDJ(0:2), DTDK(0:2)
      DOUBLE PRECISION ANGLE, DI(0:2), DJ(0:2), DK(0:2), DL(0:2)
      DOUBLE PRECISION DPDI(0:2), DPDJ(0:2), DPDK(0:2), DPDL(0:2)
      DOUBLE PRECISION FNANGC, SNANGC, SNANGCR, NANGCR
      DOUBLE PRECISION INDLDH(NLDH, 9)
      INTEGER CURPHI
      INTEGER STARTLINDH
      DOUBLE PRECISION PHI2, DPDI2(0:2), DPDJ2(0:2), DPDK2(0:2), DPDL2(0:2)
      DOUBLE PRECISION TEST(3*NATOMS,3*NATOMS)
 
      DOUBLE PRECISION PHI1, PREVPHI1

      INTEGER NBOND

      TEST(:,:) = 0.0

      NBOND = NBONA +NBONH
      IF (.NOT.ALLOCATED(IBIB).OR..NOT.ALLOCATED(JBJB)) THEN
      PRINT*, "AMBERT: BONDS NOT ALLOCATED. CHECK AMBGETNATINTERN. STOP"
      STOP
      ENDIF

      IF (NCART.NE.3*NATOMS) THEN
         PRINT*, 'GETBEENAT>> ERROR! NCART NOT 3*NATOMS', NCART, NATOMS
         STOP
      ENDIF
!---------------------------------------------------

!ZERO THE INTERNAL COORD VECTOR IF NOT NOCOORD
       IF (.NOT.NOCOOR) XINT(1:NINTC) = 0.0D0

!ZERO THE B MATRIX VALUES
       IF (.NOT.NODERV) THEN
          VAL(1:NNZ) = 0.0D0
          INDX(1:NNZ) = 0
          JNDX(1:NNZ) = 0
       ENDIF
!PUT CARTESIANS FROM XCART INTO X,Y,Z
!
      DO I1=1,NATOMS
        X(I1)=XCART(3*(I1-1)+1)
        Y(I1)=XCART(3*(I1-1)+2)
        Z(I1)=XCART(3*(I1-1)+3)
      ENDDO
!
!KD1 IS DIMENSION OF GCS SPARSE GC MATRIX
!IN FACT USING MATRIX DOUBLE THIS SIZE SO CAN CALCULATE UPPER AND LOWER TRIANGLE SEPARATELY
!WILL COMBINE THE TWO AT THE END

      KD1=KD+1
!
!ZERO GCS HERE
!
      IF (.NOT.NODERV) GCS(:,:) = 0.0D0

!
!THE FOLLOWING LOOPS WILL PUT INTERNALS IN XINT AND PASS THEM BACK
!IF THE FLAG NOCOOR IS FALSE
!AND PUT DERIVATIVES IN VAL IF NODERV IS FALSE

      CURBIND = 1 ! CURRENT INDEX IN B-MATRIX LISTS
      DIHC = 0

      MM = 0
      DO 10 TCOUNT=1,NBOND                 
        I=IBIB(TCOUNT)
        J=JBJB(TCOUNT)
        ! IGNORE BOND IF USING CARTESIANS FOR BOTH ATOMS
        IF ((USECART(I).AND.USECART(J)) ) CYCLE

        MM = MM + 1
        II=3*(I-1)+1
        II1=3*(I-1)+2
        II2=3*(I-1)+3
        JJ=3*(J-1)+1
        JJ1=3*(J-1)+2
        JJ2=3*(J-1)+3
 
!MSB50 - NO IDEA WHAT THE NEXT THREE LINES WOULD CORRESPOND TO IN AMBER, OR WHY THEY ARE THERE... PLEASE FIX THIS IF YOU CAN.
!      IC=ICB(TCOUNT)
!      IF(IC.EQ.0) GOTO 10
!      IF(CBC(IC).EQ.0) GOTO 10

        RX=X(I)-X(J)
        RY=Y(I)-Y(J)
        RZ=Z(I)-Z(J)
        S2=RX*RX + RY*RY + RZ*RZ
        S=SQRT(S2)
        IF (.NOT.NOCOOR) THEN
!DAE PUT S INTO XINT: IT IS INTERNAL COORDINATE
          XINT(MM)=S
        ENDIF
        IF (PRINTCOORDS) PRINT '(A14,I4,2(I4,A4),F7.3)', 'GETBEE> BOND', MM, I, IH(M04+I-1),J,IH(M04+J-1), XINT(MM)

        IF (.NOT.NODERV) THEN
!DAE CALCULATE B-MATRIX ELEMENTS FOR THIS BOND
!DAE VAL IS 1-D ARRAY OF NON-ZERO B ELEMENTS
!INDX IS ARRAY OF THE IC(ROW) INDICES OF VAL
!JNDX IS ARRAY OF THE CART(COLUMN) INDICES OF VAL
        SR=1/S
        RXSR=RX*SR
        RYSR=RY*SR
        RZSR=RZ*SR

        VAL(CURBIND)=RXSR
        INDX(CURBIND)=MM
        JNDX(CURBIND)=3*(I-1)+1

        TEST(MM, 3*(I-1)+1) = RXSR
!
        VAL(CURBIND+1)=RYSR
        INDX(CURBIND+1)=MM
        JNDX(CURBIND+1)=3*(I-1)+2

        TEST(MM,3*(I-1)+2) = RYSR
!
        VAL(CURBIND+2)=RZSR
        INDX(CURBIND+2)=MM
        JNDX(CURBIND+2)=3*(I-1)+3

        TEST(MM,3*(I-1)+3) = RZSR
!
        VAL(CURBIND+3)=-RXSR
        INDX(CURBIND+3)=MM
        JNDX(CURBIND+3)=3*(J-1)+1

        TEST(MM,3*(J-1)+1) = -RXSR
!
        VAL(CURBIND+4)=-RYSR
        INDX(CURBIND+4)=MM
        JNDX(CURBIND+4)=3*(J-1)+2

        TEST(MM,3*(J-1)+2) = -RYSR
!
        VAL(CURBIND+5)=-RZSR
        INDX(CURBIND+5)=MM
        JNDX(CURBIND+5)=3*(J-1)+3

        TEST(MM,3*(J-1)+3) = -RZSR
!
!UPDATE GEE FOR THIS IC
        IIKD=II+KD1
        IIKD1=II1+KD1
        IIKD2=II2+KD1
        JJKD=JJ+KD1
        JJKD1=JJ1+KD1
        JJKD2=JJ2+KD1

        VXI=VAL(CURBIND)
        VYI=VAL(CURBIND+1)
        VZI=VAL(CURBIND+2)
        GCS(KD1,II)=GCS(KD1,II)+VXI**2
        GCS(IIKD-II1,II1)=GCS(IIKD-II1,II1)+VXI*VYI
        GCS(IIKD-II2,II2)=GCS(IIKD-II2,II2)+VXI*VZI
        GCS(KD1,II1)=GCS(KD1,II1)+VYI**2
        GCS(IIKD1-II2,II2)=GCS(IIKD1-II2,II2)+VYI*VZI
        GCS(KD1,II2)=GCS(KD1,II2)+VZI**2
!
        VXJ=VAL(CURBIND+3)
        VYJ=VAL(CURBIND+4)
        VZJ=VAL(CURBIND+5)
        GCS(KD1,JJ)=GCS(KD1,JJ)+VXJ**2
        GCS(JJKD-JJ1,JJ1)=GCS(JJKD-JJ1,JJ1)+VXJ*VYJ
        GCS(JJKD-JJ2,JJ2)=GCS(JJKD-JJ2,JJ2)+VXJ*VZJ
        GCS(KD1,JJ1)=GCS(KD1,JJ1)+VYJ**2
        GCS(JJKD1-JJ2,JJ2)=GCS(JJKD1-JJ2,JJ2)+VYJ*VZJ
        GCS(KD1,JJ2)=GCS(KD1,JJ2)+VZJ**2
!
        GCS(IIKD-JJ,JJ)=GCS(IIKD-JJ,JJ)+VXI*VXJ
        GCS(IIKD-JJ1,JJ1)=GCS(IIKD-JJ1,JJ1)+VXI*VYJ
        GCS(IIKD-JJ2,JJ2)=GCS(IIKD-JJ2,JJ2)+VXI*VZJ
        GCS(IIKD1-JJ,JJ)=GCS(IIKD1-JJ,JJ)+VYI*VXJ
        GCS(IIKD1-JJ1,JJ1)=GCS(IIKD1-JJ1,JJ1)+VYI*VYJ
        GCS(IIKD1-JJ2,JJ2)=GCS(IIKD1-JJ2,JJ2)+VYI*VZJ
        GCS(IIKD2-JJ,JJ)=GCS(IIKD2-JJ,JJ)+VZI*VXJ
        GCS(IIKD2-JJ1,JJ1)=GCS(IIKD2-JJ1,JJ1)+VZI*VYJ
        GCS(IIKD2-JJ2,JJ2)=GCS(IIKD2-JJ2,JJ2)+VZI*VZJ
!
        CURBIND = CURBIND + 6
        ENDIF
   10  CONTINUE

      CURCOORD = NBDS+1
!     PRINT*, 'AFTER BONDS, CURBIND=', CURBIND

      DO 20 CNT=1,NCNT

        INFO(0:5) = CENTERS(CNT,0:5)
        CNTTYPE = INFO(0)

        IF (CNTTYPE.LE.3) THEN
           CNTATMS = 4 ! NUMBER OF ATOMS BOUND TO CENTER
           NANGC = 5 ! NUMBER OF ANGLE-SUM NATURAL COORDS
        ELSE IF (CNTTYPE.EQ.6) THEN
           CNTATMS = 2
           NANGC = 1
        ELSE IF (CNTTYPE.EQ.7) THEN
           CNTATMS = 3
           NANGC = 1
       ! ELSE IF (CNTTYPE.EQ.8) THEN
       !    CNTATMS = 4
       !    NANGC = 4
       !    CNTTYPE = 3
       ! ELSE IF (CNTTYPE.EQ.9) THEN
       !    CNTATMS = 4
       !    NANGC = 4
       !    CNTTYPE = 2
        ELSE IF (CNTTYPE.EQ.8) THEN
           CNTATMS = 4 
           NANGC =4 
        ELSE IF (CNTTYPE.EQ.4.OR.CNTTYPE.EQ.5) THEN
           CNTATMS = 3
           NANGC = 2
        ENDIF

!GET X COORDINATE INDICES IN XCART FOR EACH ATOM INVOLVED
        DO I1=1,CNTATMS+1
           XIND(I1) = 3*(INFO(I1) - 1) + 1
        ENDDO

!DAE CALCULATE B-MATRIX ELEMENTS FOR THESE ATOMS (WILSON DECIUS + CROSS CH.4)
!EFK: ADD TOGETHER THE B-MATRIX ELEMENTS FROM THE RELEVANT PRIMITIVE INTERNALS        
!GET ANGLES AND DERIVATIVES FOR ALL ANGLES WITH THE GIVEN CENTRAL ATOM
        TCOUNT=1
        DO I1=2,CNTATMS
           DO 110 K1=I1+1,CNTATMS+1
              ! IF CENTER IS PART OF RING, SKIP ANGLE WITHIN RING
              IF((CNTTYPE.EQ.7.OR.CNTTYPE.EQ.8).AND.I1.EQ.2.AND.K1.EQ.3) GOTO 110
!
              CALL AMBGETANGLE(XCART(XIND(I1):XIND(I1)+2),                 &
     &             XCART(XIND(1):XIND(1)+2), XCART(XIND(K1):XIND(K1)+2), &
     &             THETA, DTDI, DTDJ, DTDK, NOCOOR,NODERV)              
              IF (PRINTCOORDS) PRINT '(A15,3(I4,A4),F7.2)','GETBEE> &
          &ANGLE', INFO(I1), IH(M04+ INFO(I1)-1), INFO(1), IH(M04+INFO(1)-1),&
     &      INFO(K1),IH(M04+INFO(K1)-1), THETA

               DO 100 C1=1,NANGC ! COUNTS COORDINATES FOR THIS CENTER
                 CURCOEFF = COEFF(CNTTYPE, C1, TCOUNT)

                 IF (CURCOEFF.EQ.0) GOTO 100

                 ! UPDATE INTERNAL COORDS LIST
                 IF (.NOT.NOCOOR) THEN
                    XINT(CURCOORD+C1-1) = XINT(CURCOORD+C1-1)           & 
     &                   + CURCOEFF * THETA
                 ENDIF

                 IF (.NOT.NODERV) THEN
                    ! UPDATE DERIVATIVES IN B MATRIX
                    DO X1=0,2   ! COUNTS X,Y,Z COORDINATES
                       VAL(CURBIND + ATOMXPOSINC(CNTTYPE,C1,I1) + X1) = &
     &                      VAL(CURBIND + ATOMXPOSINC(CNTTYPE,C1,I1)+X1)&
     &                      + CURCOEFF * DTDI(X1)
                       VAL(CURBIND + ATOMXPOSINC(CNTTYPE,C1,1) + X1) =  &
     &                      VAL(CURBIND + ATOMXPOSINC(CNTTYPE,C1,1)+X1) &
     &                      + CURCOEFF * DTDJ(X1)
                       VAL(CURBIND + ATOMXPOSINC(CNTTYPE,C1,K1) + X1) = &
     &                      VAL(CURBIND + ATOMXPOSINC(CNTTYPE,C1,K1)+X1)&
     &                      + CURCOEFF * DTDK(X1)

                    ENDDO
                 ENDIF
 100          CONTINUE
              TCOUNT = TCOUNT + 1
 110       CONTINUE
        ENDDO

        IF (.NOT.NODERV) THEN
!    SET B MATRIX ROW AND COLUMN VALUES
           DO C1=1,NANGC
              DO A1 = 1,CNTATMS+1
                 ATOMXSHIFT = ATOMXPOSINC(CNTTYPE,C1,A1)
                 IF (.NOT.(ATOMXSHIFT.LT.0)) THEN                 
                    DO X1=0,2
                       INDX(CURBIND + ATOMXSHIFT+X1) = CURCOORD+C1-1
                       JNDX(CURBIND + ATOMXSHIFT+X1) = XIND(A1)+X1
                    ENDDO
                 ENDIF
              ENDDO
           ENDDO
       
 
!    UPDATE GEE MATRIX ELEMENTS
           DO I1=1,CNTATMS+1
              DO J1=I1,CNTATMS+1
                 DO C1=1,NANGC
                    IF(ATOMXPOSINC(CNTTYPE,C1,I1).GE.0 .AND.            & 
     &                   ATOMXPOSINC(CNTTYPE,C1,J1).GE.0) THEN
                                ! BOTH ATOMS ARE INVOLVED IN THIS COORDINATE
                       IX = XIND(I1) ! INDEX OF X COORDINATE OF ATOM I
                       JX = XIND(J1)
                                ! FIND HOW FAR OFFSET FROM CURRENT B INDEX ARE
                                ! THE DERIVATIVES WRT EACH ATOM
                       ISHIFT = ATOMXPOSINC(CNTTYPE,C1,I1) 
                       JSHIFT = ATOMXPOSINC(CNTTYPE,C1,J1)
                       
                       IF (I1.EQ.J1) THEN
                          DO X1=0,2
                             DO X2=X1,2
                                GCS(KD1+IX+X1-(JX+X2), JX+X2) =         &
     &                               GCS(KD1+IX+X1-(JX+X2),JX+X2) +     &
     &                               VAL(CURBIND+ISHIFT+X1)*VAL(CURBIND+JSHIFT+X2)
                             ENDDO
                          ENDDO
                       ELSE
                          DO X1=0,2
                             DO X2=0,2
                                GCS(KD1+IX+X1-(JX+X2), JX+X2) =         &
     &                               GCS(KD1+IX+X1-(JX+X2),JX+X2) +     &
     &                               VAL(CURBIND+ISHIFT+X1)*VAL(CURBIND+JSHIFT+X2)
                             ENDDO
                          ENDDO
                       ENDIF
                    ENDIF
                 ENDDO
              ENDDO
           ENDDO

        IF(CNTTYPE.NE.3) THEN
           CURBIND = CURBIND + ATOMXPOSINC(CNTTYPE, NANGC, CNTATMS+1)+3
        ELSE
           CURBIND = CURBIND + ATOMXPOSINC(CNTTYPE, NANGC, CNTATMS)+3
        ENDIF
      ENDIF

        CURCOORD = CURCOORD + NANGC
 20   CONTINUE

!    PRINT*, 'AFTER CNT CURBIND=', CURBIND

!!
!DEAL WITH RINGS
!!
      DO 25 RNG=1,NRNG

         IF(RINGS(RNG,6).EQ.-1) THEN
            RATMS = 5 ! 5 ATOMS IN RING
            NANGC = 2; NDIHC = 2 ! 2 ANGLE & 2 DIHEDRAL COORDS
            RTYPE = 25
         ELSE 
            RATMS = 6 ! 6 ATOMS IN RING
            NANGC = 3; NDIHC = 3            
            RTYPE = 26
         ENDIF
         
         DO I1=1,RATMS
            XIND(I1) = 3*(RINGS(RNG,I1) - 1) + 1
         ENDDO
        
         TCOUNT = 1
         DO I1 = 1,RATMS
            J1 = MOD(I1,RATMS)+1
            K1 = MOD(I1+1,RATMS)+1
            L1 = MOD(I1+2,RATMS)+1

            CALL AMBGETANGLE(XCART(XIND(I1):XIND(I1)+2),                   &
     &           XCART(XIND(J1):XIND(J1)+2), XCART(XIND(K1):XIND(K1)+2),&
     &           THETA, DTDI, DTDJ, DTDK, NOCOOR,NODERV)              

            CALL AMBGETTORSION(XCART(XIND(I1):XIND(I1)+2),                 &
     &           XCART(XIND(J1):XIND(J1)+2), XCART(XIND(K1):XIND(K1)+2),&
     &           XCART(XIND(L1):XIND(L1)+2), PHI, DPDI, DPDJ, DPDK, DPDL,&
     &           NOCOOR,NODERV)


            DIHC = DIHC + 1

            IF (TCOUNT.EQ.1) PREVPHI1 = PREVDIH(DIHC)
          
 
            IF (.NOT.NOCOOR) CALL AMBALIGNDIH(PHI, DIHC, PHI-PHI1, PREVDIH(DIHC)-PREVPHI1,TCOUNT)
            IF (TCOUNT.EQ.1) PHI1 = PHI

            DO C1 = 1,NANGC+NDIHC
               IF (C1.LE.NANGC) THEN
                  ANGLE = THETA; DI(:) = DTDI(:)
                  DJ(:)=DTDJ(:); DK(:)=DTDK(:)
               ELSE
                  ANGLE = PHI
                  DI(:)=DPDI(:); DJ(:)=DPDJ(:)
                  DK(:)=DPDK(:); DL(:)=DPDL(:)
               ENDIF

               CURCOEFF = COEFF(RTYPE, C1, TCOUNT)
               
               IF(.NOT.NOCOOR) XINT(CURCOORD+C1-1) = &
     &              XINT(CURCOORD+C1-1)+CURCOEFF*ANGLE

               IF (.NOT.NODERV) THEN
!    UPDATE B-MATRIX ELEMENTS
               DO X1=0,2        ! COUNTS X,Y,Z COORDINATES
                  VAL(CURBIND + RATMS*3*(C1-1) + 3*(I1-1) + X1) =&
     &                 VAL(CURBIND + RATMS*3*(C1-1)+ 3*(I1-1) + X1)&
     &                 + CURCOEFF * DI(X1)
                  VAL(CURBIND + RATMS*3*(C1-1) + 3*(J1-1) + X1) =&
     &                 VAL(CURBIND + RATMS*3*(C1-1) + 3*(J1-1) + X1)&
     &                 + CURCOEFF * DJ(X1)
                  VAL(CURBIND + RATMS*3*(C1-1) + 3*(K1-1) + X1) =&
     &                 VAL(CURBIND + RATMS*3*(C1-1)+ 3*(K1-1) + X1)&
     &                 + CURCOEFF * DK(X1)                  
               ENDDO
               IF (C1.GT.NANGC) THEN
                  DO X1=0,2
                     VAL(CURBIND + RATMS*3*(C1-1) + 3*(L1-1) + X1) =&
     &                    VAL(CURBIND + RATMS*3*(C1-1)+ 3*(L1-1) + X1)&
     &                    + CURCOEFF * DL(X1)
                  ENDDO
               ENDIF
               ENDIF
            ENDDO
         TCOUNT = TCOUNT + 1
         ENDDO

!SET B MATRIX ROW AND COLUMN VALUES
         IF (.NOT.NODERV) THEN
         DO C1=0,NANGC+NDIHC-1  ! COUNT COORDS
            DO A1 = 0,RATMS-1   ! COUNT ATOMS
               DO X1=0,2        ! COUNT X,Y,Z
                  INDX(CURBIND+RATMS*3*C1+3*A1+X1) = CURCOORD+C1
                  JNDX(CURBIND+RATMS*3*C1+3*A1+X1) = XIND(A1+1)+X1
               ENDDO
            ENDDO
         ENDDO
         
!UPDATE G MATRIX ELEMENTS
         DO I1=1,RATMS
            DO J1=I1,RATMS
               DO C1=1,NANGC+NDIHC

                  IX = XIND(I1)
                  JX = XIND(J1)
                  ISHIFT = RATMS*3*(C1-1)+3*(I1-1)
                  JSHIFT = RATMS*3*(C1-1)+3*(J1-1)

                  IF (I1.EQ.J1) THEN
                     DO X1=0,2
                        DO X2 = X1,2
                           GCS(KD1+IX+X1-(JX+X2), JX+X2) = &
     &                          GCS(KD1+IX+X1-(JX+X2), JX+X2) &
     &                          + VAL(CURBIND+ISHIFT+X1) * VAL(CURBIND+JSHIFT+X2)
                        ENDDO
                     ENDDO
                  ELSE
                     DO X1=0,2
                        DO X2 = 0,2
                           GCS(KD1+IX+X1-(JX+X2), JX+X2) = &
     &                          GCS(KD1+IX+X1-(JX+X2), JX+X2) &
     &                          + VAL(CURBIND+ISHIFT+X1) * VAL(CURBIND+JSHIFT+X2)                           
                        ENDDO
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
         CURBIND = CURBIND + RATMS*3*(NANGC + NDIHC)
         ENDIF
         CURCOORD = CURCOORD + NANGC + NDIHC

 25   CONTINUE
!     PRINT*, 'AFTER RNG, CURBIND:', CURBIND

!
!DEAL WITH FUSED RINGS
!GET DIFFERENCE OF DIHEDRALS 3-1-2-6 & 4-2-1-5
!
      DO FRG = 1,NFRG
         DO I1=1,6
            XIND(I1) = 3*(FRINGS(FRG,I1) - 1)+1
         ENDDO

         CALL AMBGETTORSION(XCART(XIND(3):XIND(3)+2),&
     &           XCART(XIND(1):XIND(1)+2), XCART(XIND(2):XIND(2)+2), &
     &           XCART(XIND(6):XIND(6)+2), PHI, DPDI, DPDJ, DPDK, DPDL, &
     &           NOCOOR,NODERV)   
         DIHC = DIHC + 1
         PREVPHI1 = PREVDIH(DIHC)
         IF (.NOT.NOCOOR) CALL AMBALIGNDIH(PHI, DIHC, 0.0D0, 0.0D0, 1)         

         CALL AMBGETTORSION(XCART(XIND(5):XIND(5)+2),&
     &           XCART(XIND(1):XIND(1)+2), XCART(XIND(2):XIND(2)+2),   &
     &           XCART(XIND(4):XIND(4)+2), PHI2, DPDI2, DPDJ2, DPDK2, DPDL2,&
     &           NOCOOR,NODERV)
         DIHC = DIHC + 1
         
         IF (.NOT.NOCOOR) CALL AMBALIGNDIH(PHI2, DIHC, PHI2-PHI, PREVDIH(DIHC)-PREVPHI1, 2)         
         
         IF (.NOT.NOCOOR) XINT(CURCOORD) = XINT(CURCOORD) + (PHI - PHI2)/SQRT(2.0D0)


         IF (.NOT.NODERV) THEN

!UPDATE B MATRIX
         DO X1 = 0,2
            VAL(CURBIND+X1) = VAL(CURBIND+X1) + (DPDJ(X1) - DPDJ2(X1))/SQRT(2.0D0)
            VAL(CURBIND+X1+3) = VAL(CURBIND+X1+3) + (DPDK(X1)-DPDK2(X1))/SQRT(2.0D0)
            VAL(CURBIND+X1+6) = VAL(CURBIND+X1+6) + DPDI(X1)/SQRT(2.0D0)
            VAL(CURBIND+X1+9) = VAL(CURBIND+X1+9) - DPDL2(X1)/SQRT(2.0D0)
            VAL(CURBIND+X1+12) = VAL(CURBIND+X1+12) - DPDI2(X1)/SQRT(2.0D0)
            VAL(CURBIND+X1+15) = VAL(CURBIND+X1+15) + DPDL(X1)/SQRT(2.0D0)
            DO A1=1,6
               INDX(CURBIND+3*(A1-1)+X1) = CURCOORD
               JNDX(CURBIND+3*(A1-1)+X1) = XIND(A1)+X1
            ENDDO
         ENDDO
!UPDATE GCS
         DO I1=1,6
            DO J1=I1,6
               IX = XIND(I1)
               JX = XIND(J1)

               ISHIFT = 3*(I1-1)
               JSHIFT = 3*(J1-1)

               IF (I1.EQ.J1) THEN
                  DO X1=0,2
                     DO X2=X1,2
                        GCS(KD1+IX+X1-(JX+X2), JX+X2) =        &
     &                       GCS(KD1+IX+X1-(JX+X2), JX+X2) +   &
     &                       VAL(CURBIND+ISHIFT+X1)*VAL(CURBIND+JSHIFT+X2)
                     ENDDO
                  ENDDO
               ELSE
                  DO X1=0,2
                     DO X2=0,2
                        GCS(KD1+IX+X1-(JX+X2), JX+X2) =&
     &                       GCS(KD1+IX+X1-(JX+X2), JX+X2) +&
     &                       VAL(CURBIND+ISHIFT+X1)*VAL(CURBIND+JSHIFT+X2)
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
         CURBIND = CURBIND + 18
         ENDIF
         CURCOORD = CURCOORD + 1
      ENDDO
!     PRINT*, 'AFTER FRG CURBIND:', CURBIND

!
!DEAL WITH LINEAR DIHEDRALS
!
      STARTLINDH = DIHC +1 !LINE 1272 --> FOR AMBALIGN-CALL
      DO 30 LDH=1,NLDH
         INFO(-1:8) = LINDIH(LDH,-1:8)
         !INFO(-1) AND INFO(0) DENOTE THE NUMBER OF POSSIBLE I AND L ATOMS
         DO I1=1,INFO(0)+INFO(-1)+2
            XIND(I1) = 3*(INFO(I1) - 1) + 1
         ENDDO
         NANGC = INFO(-1)*INFO(0) ! NUMBER OF TORSIONS INVOLVED

         NANGCR = 1.0D0/NANGC

         CURPHI = 0         
         DO I1=3,INFO(-1)+2
            DO L1=INFO(-1)+3,INFO(-1)+INFO(0)+2
               CURPHI = CURPHI + 1
               
               IX = XIND(I1)
               JX = XIND(1)
               KX = XIND(2)
               LX = XIND(L1)

               CALL AMBGETTORSION(XCART(IX:IX+2), XCART(JX:JX+2),  &
     &              XCART(KX:KX+2), XCART(LX:LX+2), PHI, DPDI, DPDJ,   &
     &              DPDK, DPDL, NOCOOR,NODERV)                              

               DIHC = DIHC + 1
               IF (CURPHI.EQ.1) PREVPHI1 = PREVDIH(DIHC)
 
               IF (PRINTCOORDS) PRINT '(A15,3I4,4(I4,A4),F8.2)','GETBEE>& 
            &DIH',LDH, CURCOORD,DIHC, INFO(I1), IH(M04+INFO(I1)-1),INFO(1),IH(M04+INFO(1)-1),&
              INFO(2), IH(M04+INFO(2)-1),INFO(L1),IH(M04+INFO(L1)-1), PHI
               
               IF (.NOT.NOCOOR) CALL AMBALIGNDIH(PHI, DIHC, PHI-PHI1, PREVDIH(DIHC)-PREVPHI1,CURPHI)

               IF (CURPHI.EQ.1) PHI1 = PHI                  

               IF (.NOT.NOCOOR) XINT(CURCOORD) = XINT(CURCOORD) + PHI*NANGCR
               IF (.NOT.NODERV) THEN
               ! UPDATE B-MATRIX ELEMENTS WITH THIS TORSION
               DO X1=0,2
                  VAL(CURBIND+3*(I1-1)+X1)=VAL(CURBIND+3*(I1-1)+X1) + DPDI(X1)*NANGCR
                  VAL(CURBIND+X1) = VAL(CURBIND+X1) + DPDJ(X1)*NANGCR
                  VAL(CURBIND+X1+3) = VAL(CURBIND+X1+3) + DPDK(X1)*NANGCR
                  VAL(CURBIND+3*(L1-1)+X1)=VAL(CURBIND+3*(L1-1)+X1) + DPDL(X1)*NANGCR
               ENDDO
               ENDIF
            ENDDO
         ENDDO
          
         IF (.NOT.NODERV) THEN
         ! SET B-MATRIX ROW AND COLUMN INDICES
         DO A1=1,INFO(0)+INFO(-1)+2
            DO X1=0,2
               INDX(CURBIND+3*(A1-1)+X1) = CURCOORD
               JNDX(CURBIND+3*(A1-1)+X1) = XIND(A1) + X1
            ENDDO
         ENDDO

         ! UPDATE THE GEE MATRIX ELEMENTS
         DO I1=1,INFO(0)+INFO(-1)+2
            DO J1=I1,INFO(0)+INFO(-1)+2
               IX = XIND(I1)
               JX = XIND(J1)
               IF(I1.EQ.J1) THEN
                  DO X1=0,2
                     DO X2=X1,2
                        GCS(KD1+IX+X1-(JX+X2), JX+X2) = &
     &                       GCS(KD1+IX+X1-(JX+X2), JX+X2) &
     &                      + VAL(CURBIND+3*(I1-1)+X1)*VAL(CURBIND+3*(J1-1)+X2)
                     ENDDO
                  ENDDO
               ELSE
                  DO X1=0,2
                     DO X2=0,2
                        GCS(KD1+IX+X1-(JX+X2), JX+X2) =&
     &                       GCS(KD1+IX+X1-(JX+X2), JX+X2)&
     &                      + VAL(CURBIND+3*(I1-1)+X1)*VAL(CURBIND+3*(J1-1)+X2)
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
         CURBIND = CURBIND + 3*(INFO(0) + INFO(-1) + 2)
         ENDIF
         CURCOORD = CURCOORD + 1
         

   30 CONTINUE
      !PRINT*, 'AFTER LDH CURBIND=', CURBIND

      DO 40 IMP=1,NIMP
!IMPROPER DIHEDRAL ('WAGGING') COORDINATES
         INFO(1:4) = IMPDIH(IMP,1:4)
         IX = 3*(INFO(1)-1)+1
         JX = 3*(INFO(2)-1)+1
         KX = 3*(INFO(3)-1)+1
         LX = 3*(INFO(4)-1)+1

         CALL AMBGETOUTOFPLANE(XCART(IX:IX+2), XCART(JX:JX+2),&
     &              XCART(KX:KX+2), XCART(LX:LX+2), PHI, DPDI, DPDJ,&
     &              DPDK, DPDL, NOCOOR,NODERV)


         IF(.NOT.NOCOOR) XINT(CURCOORD) = PHI

         IF (.NOT.NODERV) THEN
!GET B-MATRIX ELEMENTS
         DO X1=0,2
            VAL(CURBIND+X1) = DPDI(X1)
            INDX(CURBIND+X1) = CURCOORD
            JNDX(CURBIND+X1) = IX+X1

            VAL(CURBIND+X1+3) = DPDJ(X1)
            INDX(CURBIND+X1+3) = CURCOORD
            JNDX(CURBIND+X1+3) = JX+X1
            
            VAL(CURBIND+X1+6) = DPDK(X1)
            INDX(CURBIND+X1+6) = CURCOORD
            JNDX(CURBIND+X1+6) = KX+X1

            VAL(CURBIND+X1+9) = DPDL(X1)
            INDX(CURBIND+X1+9) = CURCOORD
            JNDX(CURBIND+X1+9) = LX+X1

         ENDDO

!UPDATE GEE MATRIX
         DO I1=1,4
            DO J1=I1,4
               IX = 3*(INFO(I1)-1)+1
               JX = 3*(INFO(J1)-1)+1
               IF (I1.EQ.J1) THEN
                  DO X1=0,2
                     DO X2=X1,2                     
                        GCS(KD1+IX+X1-(JX+X2), JX+X2) =                 &    
     &                       GCS(KD1+IX+X1-(JX+X2), JX+X2)              &
     &                       + VAL(CURBIND+3*(I1-1)+X1)*VAL(CURBIND+3*(J1-1)+X2)
                     ENDDO
                  ENDDO
               ELSE
                  DO X1=0,2
                     DO X2=0,2                     
                        GCS(KD1+IX+X1-(JX+X2), JX+X2) =                 &
     &                       GCS(KD1+IX+X1-(JX+X2), JX+X2)              &
     &                       + VAL(CURBIND+3*(I1-1)+X1)*VAL(CURBIND+3*(J1-1)+X2)                        
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
         CURBIND = CURBIND + 12
         ENDIF
         CURCOORD = CURCOORD + 1
         
   40 CONTINUE
!     PRINT*, 'AFTER IMP CURBIND=', CURBIND

      DO 50 MM = 1,NCRT
         I = CARTATMS(MM)
         IF (.NOT.NOCOOR) XINT(CURCOORD:CURCOORD+2) = XCART(3*(I-1)+1:3*I)

         IF (.NOT.NODERV) THEN
            DO X1 = 0,2
               VAL(CURBIND+X1) = 1.0D0
               INDX(CURBIND+X1) = CURCOORD+X1
               JNDX(CURBIND+X1) = 3*(I-1)+1+X1
               GCS(KD1,3*(I-1)+1+X1) = GCS(KD1,3*(I-1)+1+X1) + 1.0D0
            ENDDO
            CURBIND = CURBIND + 3
         ENDIF

         CURCOORD = CURCOORD + 3         
 50   CONTINUE
!     PRINT*, 'AFTER CART CURBIND=', CURBIND

      IF (.NOT.NODERV) THEN
!ADD LOWER DIAGONAL ELEMENTS TO UPPER DIAGONAL. AVOID DOUBLING DIAGONAL
!BY LOOPING ONLY TO J1-1
       DO J1=1,KD
         DO I1=1,J1-1
           GCS(KD+1+I1-J1,J1)=GCS(KD+1+I1-J1,J1)+GCS(KD+1+J1-I1,I1)
         ENDDO
       ENDDO
       DO J1=KD+1,NCART
         DO I1=J1-KD,J1-1
           GCS(KD+1+I1-J1,J1)=GCS(KD+1+I1-J1,J1)+GCS(KD+1+J1-I1,I1)
         ENDDO
       ENDDO
      ENDIF
      
      IF (PRINTCOORDS) THEN
         PRINT*, 'DONE PRINTING COORDS'
         PRINT*, XINT
!         STOP
      ENDIF

      RETURN
      END SUBROUTINE


