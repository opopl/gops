!   CONNECT MODULE IS AN IMPLEMENTATION OF A CONNECTION ALGORITHM FOR FINDING REARRANGEMENT PATHWAYS.
!   COPYRIGHT (C) 2003-2006 SEMEN A. TRYGUBENKO AND DAVID J. WALES
!   THIS FILE IS PART OF CONNECT MODULE. CONNECT MODULE IS PART OF OPTIM.
!
!   OPTIM IS FREE SOFTWARE; YOU CAN REDISTRIBUTE IT AND/OR MODIFY
!   IT UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE AS PUBLISHED BY
!   THE FREE SOFTWARE FOUNDATION; EITHER VERSION 2 OF THE LICENSE, OR
!   (AT YOUR OPTION) ANY LATER VERSION.
!
!   OPTIM IS DISTRIBUTED IN THE HOPE THAT IT WILL BE USEFUL,
!   BUT WITHOUT ANY WARRANTY; WITHOUT EVEN THE IMPLIED WARRANTY OF
!   MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  SEE THE
!   GNU GENERAL PUBLIC LICENSE FOR MORE DETAILS.
!
!   YOU SHOULD HAVE RECEIVED A COPY OF THE GNU GENERAL PUBLIC LICENSE
!   ALONG WITH THIS PROGRAM; IF NOT, WRITE TO THE FREE SOFTWARE
!   FOUNDATION, INC., 59 TEMPLE PLACE, SUITE 330, BOSTON, MA  02111-1307  USA
!
MODULE CONNECTUTILS
     USE CONNECTDATA
     IMPLICIT NONE
     INTEGER NREDO
     CONTAINS

     FUNCTION ISNEWTS(TSTOCHECK)
          USE NEBTOCONNECT
          USE KEYUTILS
          USE KEY,ONLY : DEBUG, RIGIDBODY, BULKT, TWOD, PERMDIST
          USE COMMONS,ONLY : PARAM1,PARAM2,PARAM3
          IMPLICIT NONE
          LOGICAL :: ISNEWTS
          TYPE(TSFOUNDTYPE),INTENT(IN) :: TSTOCHECK
          DOUBLE PRECISION RMAT(3,3), DIST2
          
          INTEGER :: I

          IF (NTS==0) THEN
               ISNEWTS=.TRUE.
               RETURN
          ENDIF

          DO I=1,NTS
               IF (ABS(TSTOCHECK%E-TS(I)%DATA%E) < EDIFFTOL) THEN
!PRINT *, "ENERGY OF THE TS FOUND IS THE SAME AS OF THE TS #",I,"; CHECKING DISTANCE..."
                  CALL MINPERMDIST(TSTOCHECK%COORD,TS(I)%DATA%X, NATOMS, &
  &                                DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,D,DIST2,RIGIDBODY,RMAT)
!PRINT *, "THE DISTANCE IS",D
                    IF (D<GEOMDIFFTOL) THEN
                         !PRINT *, "THIS TS IS ALREADY KNOWN."
                         ISNEWTS=.FALSE.
!                        IF (I > NTSOLD) THEN ! TWO OR MORE TS RETURNED FROM _ONE_ NEB RUN ARE THE SAME
!                             PRINT '(A,I8)', ' ISNEWTS> SHORTCUT WAS FOUND - TS IS SAME AS ALREADY SAVED TS #',I 
!                        ENDIF
                         PRINT '(A,I8,A,F20.10)', ' ISNEWTS> TRANSITION STATE IS SAME AS NUMBER ',I,' ENERGY=',TSTOCHECK%E
                         RETURN
                    ENDIF
               ENDIF
          ENDDO

          ISNEWTS=.TRUE.
     END FUNCTION ISNEWTS
    
     SUBROUTINE ISNEWMIN(E,COORD,MINPOS,NEW,REDOPATH,PERMUTE,INVERT,INDEX,I)
          USE KEYUTILS
          USE KEY,ONLY : DEBUG, RIGIDBODY, BULKT, TWOD, PERMDIST, AMBERT !MSB50
          USE COMMONS,ONLY : PARAM1,PARAM2,PARAM3
          USE INTCOMMONS, ONLY : INTMINPERMT, INTINTERPT
          IMPLICIT NONE
          
          DOUBLE PRECISION,POINTER     :: E,COORD(:) 
          INTEGER,INTENT(OUT) :: MINPOS ! IF MINIMUM NEW IT IS EQUAL TO NMIN+1; OTHERWISE - POSITION IN MIN ARRAY
          LOGICAL,INTENT(OUT) :: NEW
          DOUBLE PRECISION QTEMP(3*NATOMS)
          INTEGER :: I
          INTEGER INVERT, INDEX(NATOMS), J2
          LOGICAL PERMUTE,SUCCESS,REDOPATH
          DOUBLE PRECISION D, DIST2, RMAT(3,3)

          MINPOS=NMIN+1
          NEW=.TRUE.

          DO I=1,NMIN
             PERMUTE=.FALSE.
!            IF (DEBUG) PRINT '(A,I6,2G20.10)',' ISNEWMIN> I,E,MI(I)%DATA%E=',I,E,MI(I)%DATA%E
             IF (ABS(E-MI(I)%DATA%E) < EDIFFTOL) THEN
                CALL MINPERMDIST(COORD,MI(I)%DATA%X(1:3*NATOMS), NATOMS, &
  &                                DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,D,DIST2,RIGIDBODY,RMAT)
!               IF (DEBUG) PRINT '(A,G20.10)','ISNEWMIN> MINIMUM DISTANCE=',D
                IF (D<GEOMDIFFTOL) THEN
                   NEW=.FALSE.
                   MINPOS=I
!                  IF (DEBUG) PRINT '(A,I6)','ISNEWMIN> MINPOS=',I
                   RETURN
                ENDIF
             ENDIF
          ENDDO

!                 SUCCESS=.FALSE. 
! !
! !  IF THEY ARE THE SAME MINIMUM THEN GETPERM SHOULD BE MORE RELIABLE THAN MINDIST. DJW
! !  EVEN IF THE PERMUTATION-INVERSION IS THE IDENTITY!
! !  GETPERM CHANGES THE FIRST ARGUMENT, BUT NOT THE SECOND.
! !  HOWEVER! GETPERM WILL FAIL IF THE STATIONARY POINTS ARE WITHIN EDIFFTOL BUT
! !  DON;T LINE UP SUFFICIENTLY WELL. 
! !
!                 QTEMP(1:3*NATOMS)=MI(I)%DATA%X(1:3*NATOMS)
!                 CALL GETPERM(QTEMP,COORD,INVERT,INDEX,SUCCESS)
! 
!                 IF (SUCCESS) THEN
!                    ATOMLOOP: DO J2=1,NATOMS
!                       IF ((INVERT.NE.1).OR.(INDEX(J2).NE.J2)) THEN
!                          PERMUTE=.TRUE.
!                          MINPOS=I
!                          EXIT ATOMLOOP
!                       ENDIF
!                    ENDDO ATOMLOOP
!                    IF (.NOT.PERMUTE) THEN
!                       MINPOS=I   ! SAME PERMUTATION-INVERSION ISOMER
!                       NEW=.FALSE.
!                       RETURN
!                    ELSE ! DON;T RETURN UNLESS PERMDIST IS TRUE - WE MAY NEED TO CHECK THE OTHER MINIMA FOR AN EXACT MATCH
!                       IF (INVERT.EQ.1) PRINT '(A,I6,A)',   &
!                              ' ISNEWMIN> PERMUTATIONAL ISOMER OF MINIMUM ',I,' IDENTIFIED'
!                       IF (INVERT.EQ.-1) PRINT '(A,I6,A)',  &
!   &                          ' ISNEWMIN> PERMUTATION-INVERSION ISOMER OF MINIMUM ',I,' IDENTIFIED'
!                       IF (PERMDIST) THEN
!                          MINPOS=I   ! WE SHOULD NOT NEED TO SAVE ANY PERMUTATION-INVERSION ISOMERS IF PERMDIST IS TRUE
!                          NEW=.FALSE.
!                          RETURN
!                       ENDIF
!                    ENDIF
!                 ENDIF
!              ENDIF
!           ENDDO
! 
!           IF ((MINPOS.GT.0).AND.REDOPATH) THEN
!              NEW=.FALSE. ! MINPOS WAS SET IN THE LAST MATCH
!              PERMUTE=.TRUE.
!           ELSE
!              MINPOS=NMIN+1
!              NEW=.TRUE.
!              PERMUTE=.FALSE.
!           ENDIF
     END SUBROUTINE ISNEWMIN

     SUBROUTINE TESTSAMEMIN(E,COORD,E2,COORD2,SAME)
          USE KEYUTILS
          USE KEY,ONLY : DEBUG, RIGIDBODY, BULKT, TWOD
          IMPLICIT NONE
          DOUBLE PRECISION RMAT(3,3)

          DOUBLE PRECISION,POINTER     :: E,COORD(:),E2,COORD2(:)
          LOGICAL,INTENT(OUT) :: SAME

          SAME = .FALSE.

          IF (ABS(E-E2) < EDIFFTOL) THEN
!                   PRINT '(A)','CALLING MINDIST 3'
              CALL NEWMINDIST(COORD,COORD2,NATOMS,D,BULKT,TWOD,'AX   ',.TRUE.,RIGIDBODY,DEBUG,RMAT)
              IF (D < GEOMDIFFTOL) THEN
                  SAME = .TRUE.
                  WRITE(*,'(A)') ' SAME MINIMUM FOUND ON EITHER SIDE OF THE PATHWAY - REJECTING THIS MIN--SAD--MIN TRIPLE'
              ENDIF
          ENDIF
     END SUBROUTINE TESTSAMEMIN

     SUBROUTINE ADDNEWMIN(E,COORD)
          USE KEYDECIDE,ONLY : INTERPCOSTFUNCTION
          USE KEY,ONLY : DEBUG, RIGIDBODY, PERMDIST,TWOD,BULKT,INTCONSTRAINTT,INTLJT,INTIMAGE,FREEZENODEST, &
  &                      ATOMACTIVE
          USE COMMONS,ONLY : PARAM1,PARAM2,PARAM3
          IMPLICIT NONE
          DOUBLE PRECISION RMAT(3,3), DIST2
          DOUBLE PRECISION,POINTER :: E,COORD(:)
          DOUBLE PRECISION VNEW(3*NATOMS), RMS, ENERGY  !!! DJW
          INTEGER :: I, NMAXINT, NMININT, INTIMAGESAVE, J1, J2
          DOUBLE PRECISION MINCOORDS(2,3*NATOMS), CONSTRAINTE, XYZLOCAL(6*NATOMS), GGGLOCAL(6*NATOMS), RMSLOCAL, &
  &                        EEELOCAL(INTIMAGE+2)
          LOGICAL IMGFREEZELOCAL(0), FREEZENODESTLOCAL
!         LOGICAL EDGEINT(INTIMAGE+1,NATOMS,NATOMS)

          IF (NMIN==MINRACKSIZE) CALL REALLOCATEMINRACK

          NMIN=NMIN+1
          MI(NMIN)%DATA%E => E
          MI(NMIN)%DATA%X => COORD
!         PRINT '(A,I6)',' ADDNEWMIN> NMIN=',NMIN
!         DO I=1,NMIN
!           PRINT '(A,I6,F20.10)',' ADDNEWMIN> I,ENERGY=',I,MI(I)%DATA%E
!         ENDDO
          ALLOCATE( MI(NMIN)%DATA%D(NMIN-1),MI(NMIN)%DATA%NTRIES(NMIN-1),MI(NMIN)%DATA%C(1) )
          IF (INTERPCOSTFUNCTION) ALLOCATE( MI(NMIN)%DATA%INTERP(NMIN-1))

          DO I=1,NMIN-1
               CALL MINPERMDIST(MI(NMIN)%DATA%X(:), MI(I)%DATA%X(:), NATOMS, &
  &                             DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,D,DIST2,RIGIDBODY,RMAT)
!              IF (DEBUG) PRINT '(A,G20.10,2(A,I6))',' ADDNEWMIN> DISTANCE FROM MINPERMDIST/NEWMINDIST=',D,' FOR ',I,' AND ',NMIN
               MI(NMIN)%DATA%D(I)=D  !  PASSING MIN(NMIN)%DATA%D(I) INSTEAD OF D DOES NOT WORK
                                     !  PERHAPS BECAUSE MIN IS AN INTRISIC FUNCTION!
               IF (INTERPCOSTFUNCTION) THEN
                  IF (INTCONSTRAINTT) THEN
                     MINCOORDS(1,1:NOPT)=MI(NMIN)%DATA%X(1:NOPT)
                     MINCOORDS(2,1:NOPT)=MI(I)%DATA%X(1:NOPT)
                     CALL MAKE_CONPOT(2,MINCOORDS)
!
! NMAXINT AND NMININT ARE RETURNED.
! 
                     FREEZENODESTLOCAL=FREEZENODEST
                     FREEZENODEST=.FALSE.
                     XYZLOCAL(1:NOPT)=MINCOORDS(1,1:NOPT)
                     XYZLOCAL(NOPT+1:2*NOPT)=MINCOORDS(2,1:NOPT)
                     INTIMAGESAVE=INTIMAGE
                     INTIMAGE=0
                     CALL CONGRAD2(NMAXINT,NMININT,CONSTRAINTE,XYZLOCAL,GGGLOCAL,EEELOCAL,IMGFREEZELOCAL,RMSLOCAL)
                     MI(NMIN)%DATA%INTERP(I)=CONSTRAINTE/2.0D0 ! ENERGY PER IMAGE
                     IF (DEBUG) PRINT '(A,I6,A,I6,2(A,G20.10))',' ADDNEWMIN> INTERPOLATION METRIC FOR ',NMIN,&
  &                  ' AND ',I,' IS ',MI(NMIN)%DATA%INTERP(I),' DIST=',MI(NMIN)%DATA%D(I)
                     INTIMAGE=INTIMAGESAVE
                     FREEZENODEST=FREEZENODESTLOCAL
                  ELSE IF (INTLJT) THEN
                     MINCOORDS(1,1:NOPT)=MI(NMIN)%DATA%X(1:NOPT)
                     MINCOORDS(2,1:NOPT)=MI(I)%DATA%X(1:NOPT)
                     FREEZENODESTLOCAL=FREEZENODEST
                     FREEZENODEST=.FALSE.
                     XYZLOCAL(1:NOPT)=MINCOORDS(1,1:NOPT)
                     XYZLOCAL(NOPT+1:2*NOPT)=MINCOORDS(2,1:NOPT)
                     INTIMAGESAVE=INTIMAGE
                     INTIMAGE=0
                     ATOMACTIVE(1:NATOMS)=.TRUE.
!                    EDGEINT(1:INTIMAGE+1,1:NATOMS,1:NATOMS)=.FALSE.
                     CALL INTGRADLJ(CONSTRAINTE,XYZLOCAL,GGGLOCAL,IMGFREEZELOCAL,RMSLOCAL,.FALSE.)
                     MI(NMIN)%DATA%INTERP(I)=CONSTRAINTE/2.0D0 ! ENERGY PER IMAGE
                     IF (DEBUG) PRINT '(A,I6,A,I6,2(A,G20.10))',' ADDNEWMIN> INTERPOLATION METRIC FOR ',NMIN,&
  &                  ' AND ',I,' IS ',MI(NMIN)%DATA%INTERP(I),' DIST=',MI(NMIN)%DATA%D(I)
                     INTIMAGE=INTIMAGESAVE
                     FREEZENODEST=FREEZENODESTLOCAL
                  ELSE
                     MI(NMIN)%DATA%INTERP(I)=INTERPVALUE(MI(NMIN)%DATA%X(:), MI(I)%DATA%X(:),D)  
                  ENDIF
               ENDIF
          ENDDO
          MI(NMIN)%DATA%NTRIES=0
          MI(NMIN)%DATA%C = 0
          MI(NMIN)%DATA%S = .FALSE.
          MI(NMIN)%DATA%F = .FALSE.
          NULLIFY( MI(NMIN)%DATA%CTS,MI(NMIN)%DATA%CMIN )
     END SUBROUTINE ADDNEWMIN

     SUBROUTINE DUMPTS
          USE NEBTOCONNECT
          USE COMMONS, ONLY: ZSYM
          IMPLICIT NONE

          INTEGER :: I,J
          OPEN(UNIT=11,FILE="TS.XYZ")
          DO I=1,NTS !FOUND 
               WRITE(11,'(I3)') NOPT/3 
               WRITE(11,'(A,I3)') 'TS #', I
               DO J=1,NOPT,3
                    WRITE(11,'(A5,1X,3F20.10)') ZSYM((J+2)/3),TS(I)%DATA%X(J),TS(I)%DATA%X(J+1),TS(I)%DATA%X(J+2)
               ENDDO
           ENDDO
           CLOSE(11)
     END SUBROUTINE DUMPTS
     
     SUBROUTINE ADDNEWCONNECTION(I,T,M) ! UPDATES CTS AND CMIN ARRAYS OF MIN I TO INCLUDE TS T AND MINIMUM M
          IMPLICIT NONE
          INTEGER,INTENT(IN) :: I,T,M
          
          INTEGER :: EXTENT
          INTEGER,POINTER :: P(:)
          
          IF ( ASSOCIATED(MI(I)%DATA%CTS) ) THEN
               EXTENT = MI(I)%DATA%C(1)
               
               !CONNECTED TS
               P => MI(I)%DATA%CTS
               NULLIFY(MI(I)%DATA%CTS)
               ALLOCATE( MI(I)%DATA%CTS(EXTENT+1) )
               MI(I)%DATA%CTS(:EXTENT) = P
               MI(I)%DATA%CTS(EXTENT+1)= T
               DEALLOCATE(P)

               !CONNECTED MINIMA
               P => MI(I)%DATA%CMIN
               NULLIFY(MI(I)%DATA%CMIN)
               ALLOCATE( MI(I)%DATA%CMIN(EXTENT+1) )
               MI(I)%DATA%CMIN(:EXTENT) = P
               MI(I)%DATA%CMIN(EXTENT+1)= M
               DEALLOCATE(P)

               MI(I)%DATA%C(1) = MI(I)%DATA%C(1) + 1
          ELSE
               ALLOCATE( MI(I)%DATA%CTS(1),MI(I)%DATA%CMIN(1) )
               MI(I)%DATA%CTS(1) = T
               MI(I)%DATA%CMIN(1)= M

               MI(I)%DATA%C = 1
          ENDIF
     END SUBROUTINE ADDNEWCONNECTION

     SUBROUTINE NEWCONNECTION(I,J,K) ! NEW CONNECTION WAS FOUND BETWEEN MINIMA I AND J (ALREADY SAVED IN MIN) VIA TS K
          USE KEYCONNECT
          IMPLICIT NONE
          INTEGER,INTENT(IN) :: I,J,K
          
          ! STORING NEW CONNECTION FOR TS AND MINIMA INVOLVED
          TS(K)%DATA%P=I
          TS(K)%DATA%M=J
          CALL ADDNEWCONNECTION(I,K,J)
          CALL ADDNEWCONNECTION(J,K,I)

          ! SETTING FLAGS:
          IF (.NOT.MI(I)%DATA%S .AND. .NOT.MI(I)%DATA%F .AND. .NOT.MI(J)%DATA%S .AND. .NOT.MI(J)%DATA%F) THEN
               PRINT *, 'CONNECTION ESTABLISHED BETWEEN MEMBERS OF THE U SET.'
               RETURN 
          ELSEIF ( (MI(I)%DATA%S.AND.MI(J)%DATA%S).OR.(MI(I)%DATA%F.AND.MI(J)%DATA%F) ) THEN
               IF (MI(I)%DATA%S) THEN
                    PRINT *, 'ALTERNATIVE PATH FOUND BETWEEN MEMBERS OF THE S SET.'
               ELSE
                    PRINT *, 'ALTERNATIVE PATH FOUND BETWEEN MEMBERS OF THE F SET.'
               ENDIF
               RETURN
          ELSEIF ( (MI(I)%DATA%S.AND.MI(J)%DATA%F).OR.(MI(J)%DATA%S.AND.MI(I)%DATA%F) ) THEN
               IF (FINISHED) RETURN ! MUST NOT CALL GETPATHWAY TWICE - SEGMENTATION FAULT! DJW
               FINISHED     =.TRUE.
               CALL GETPATHWAY(I,J)
               RETURN
          ELSEIF ( (.NOT.MI(I)%DATA%S .AND. .NOT.MI(I)%DATA%F .AND. MI(J)%DATA%S) .OR. &
                 & (.NOT.MI(J)%DATA%S .AND. .NOT.MI(J)%DATA%F .AND. MI(I)%DATA%S) ) THEN
               IF (MI(J)%DATA%S) THEN
                    WRITE(CHR,'(I7)') I
                    CALL UPDATELINK(I,J)
               ELSE
                    WRITE(CHR,'(I7)') J
                    CALL UPDATELINK(J,I)
               ENDIF
               WRITE(*,'(1X,A)') 'UNCONNECTED MINIMUM '//TRIM(ADJUSTL(CHR))//' FOUND ITS WAY TO S SET.'
          ELSEIF ( (.NOT.MI(I)%DATA%S .AND. .NOT.MI(I)%DATA%F .AND. MI(J)%DATA%F) .OR. &
                 & (.NOT.MI(J)%DATA%S .AND. .NOT.MI(J)%DATA%F .AND. MI(I)%DATA%F) ) THEN
               IF (MI(J)%DATA%F) THEN
                    WRITE(CHR,'(I7)') I
                    CALL UPDATELINK(I,J)
               ELSE
                    WRITE(CHR,'(I7)') J
                    CALL UPDATELINK(J,I)
               ENDIF
               WRITE(*,'(1X,A)') 'UNCONNECTED MINIMUM '//TRIM(ADJUSTL(CHR))//' FOUND ITS WAY TO F SET.'
          ENDIF 
     END SUBROUTINE NEWCONNECTION

     ! SEMEN Срд Лип 19 14:37:58 BST 2006: UPDATELINK IS A NEW VERSION OF UPDATELINKOLD;
     ! SHOULD BE FASTER, EASIER TO MAINTAIN AND DEBUG; UPDATELINKOLD USES THE SAME DATA STRUCTURES
     ! SO SHOULD STILL WORK --- TO ENABLE SWAP UPDATELINK AND UPDATELINKOLD.
     RECURSIVE SUBROUTINE UPDATELINK(I,J)
          USE CONNECTDATA,ONLY:LISTLENGTH
          IMPLICIT NONE
          INTEGER,INTENT(IN) :: I,J ! MINIMA I SHOULD BE BROUGHT FROM U SET TO SET TO WHICH BELONGS MINIMUM J;
          INTEGER :: K,M

          DEPTH = DEPTH+1 ! DEPTH LEVEL OF THE RECURSION
!         PRINT '(A,3I8)','UPDATELINK> I,J,DEPTH=',I,J,DEPTH
       
          ! ADD J TO THE LIST OF MINIMA TO BE AVOIDED
          IF (LISTLENGTH==0) THEN 
               ALLOCATE(START2)
               NULLIFY(START2%NEXT)
               DUMMY2 => START2
          ELSE
               ALLOCATE(DUMMY2%NEXT)
               DUMMY2 => DUMMY2%NEXT
               NULLIFY(DUMMY2%NEXT)
          ENDIF
          LISTLENGTH=LISTLENGTH+1
          DUMMY2%I=J

          ! UPDATE FLAGS
          MI(I)%DATA%S = MI(J)%DATA%S
          MI(I)%DATA%F = MI(J)%DATA%F
          !WRITE(*,'(A)',ADVANCE='NO') '.'

          ! UPDATE ALL CONNECTIONS OF MINIMUM I
       U: DO K=1,MI(I)%DATA%C(1)

               ! DO NOT TRY CONNECTIONS TO MINIMA IN THE LIST
               DUMMY2 => START2
               IF (MI(I)%DATA%CMIN(K)==DUMMY2%I) CYCLE U
               DO M=2,LISTLENGTH
                    DUMMY2 => DUMMY2%NEXT
                    IF (MI(I)%DATA%CMIN(K)==DUMMY2%I) CYCLE U
               ENDDO

               CALL UPDATELINK(MI(I)%DATA%CMIN(K),I)
          ENDDO U
          
          IF (DEPTH==1) THEN
               DO M=1,LISTLENGTH
                    DUMMY2=>START2%NEXT
                    NULLIFY(START2%NEXT)
                    DEALLOCATE(START2)
                    START2=>DUMMY2
               ENDDO
               LISTLENGTH=0
          ENDIF

          DEPTH = DEPTH-1
     END SUBROUTINE UPDATELINK

     RECURSIVE SUBROUTINE UPDATELINKOLD(I,J)
          IMPLICIT NONE
          INTEGER,INTENT(IN)          :: I,J ! MINIMA I SHOULD BE BROUGHT FROM U SET TO SET TO WHICH BELONGS MINIMUM J;
                                             ! ALL CONNECTIONS MUST BE UPDATED; CONNECTION TO MINIMUM J IS IGNORED
          INTEGER :: K,M
       
          DEPTH = DEPTH+1 ! DEPTH LEVEL OF THE RECURSION; EQUALS TO THE LENGTH OF THE LIST OF MINIMA TO BE IGNORED
!         PRINT '(A,3I8)','UPDATELINK> I,J,DEPTH=',I,J,DEPTH

          ! ADD J TO THE LIST OF MINIMA TO BE AVOIDED
          IF (DEPTH==1) THEN 
               ALLOCATE(START2)
               NULLIFY(START2%NEXT)
               DUMMY2 => START2
          ELSE
               ALLOCATE(DUMMY2%NEXT)
               DUMMY2 => DUMMY2%NEXT
               NULLIFY(DUMMY2%NEXT)
          ENDIF
          DUMMY2%I=J

          ! UPDATE FLAGS
          MI(I)%DATA%S = MI(J)%DATA%S
          MI(I)%DATA%F = MI(J)%DATA%F
          !WRITE(*,'(A)',ADVANCE='NO') '.'

          ! UPDATE ALL CONNECTIONS OF MINIMUM I
       U: DO K=1,MI(I)%DATA%C(1)

               ! DO NOT TRY CONNECTIONS TO MINIMA IN THE LIST
               DUMMY2 => START2
               IF (MI(I)%DATA%CMIN(K)==DUMMY2%I) CYCLE
               DO M=2,DEPTH
                    DUMMY2 => DUMMY2%NEXT
                    IF (MI(I)%DATA%CMIN(K)==DUMMY2%I) CYCLE U
               ENDDO

               CALL UPDATELINK(MI(I)%DATA%CMIN(K),I)
          ENDDO U
          
          ! REMOVE LAST MINIMUM FROM THE LIST
          IF (DEPTH==1) THEN
               NULLIFY(DUMMY2)
               DEALLOCATE(START2)
          ELSE
               DUMMY2=>START2
               DO K=2,DEPTH-1
                    DUMMY2 => DUMMY2%NEXT
               ENDDO
               DEALLOCATE(DUMMY2%NEXT)
          ENDIF

          DEPTH = DEPTH-1
     END SUBROUTINE UPDATELINKOLD

     SUBROUTINE GETPATHWAY(I,J)
          IMPLICIT NONE
          
          INTEGER,INTENT(IN) :: I,J ! MINIMA I AND J ARE CONNECTED TOGETHER, AND EACH WITH ONE OF THE ENDPOINTS
                                    ! YET I BELONGS TO S, J - TO F (OR VICE VERSA)
          ! CONSTRUCTING THE CHAIN
          ENDREACHED = .FALSE.
          IF (MI(I)%DATA%S) THEN
               CALL EXTRACT(1,I,J)
               ENDREACHED = .FALSE.
               CALL EXTRACT(J,2,I)
          ELSE
               CALL EXTRACT(1,J,I)
               ENDREACHED = .FALSE.
               CALL EXTRACT(I,2,J)
          ENDIF
     END SUBROUTINE GETPATHWAY

     RECURSIVE SUBROUTINE EXTRACT(H,I,J) ! I IS CONNECTED TO (THE ENDPOINT) H; SUBROUTINE EXTRACTS THIS PART OF CHAIN
          IMPLICIT NONE                  ! IN THE FORM:  START(H)--()--...-()--DUMMY(I)
          INTEGER,INTENT(IN) :: H,I,J    ! J CONNECTION IS TO BE IGNORED

          INTEGER :: K,M

          DEPTH = DEPTH+1 
          IF (DEPTH==1) THEN
               ALLOCATE(START2)
               NULLIFY(START2%NEXT)
               DUMMY2 => START2
          ELSE
               ALLOCATE(DUMMY2%NEXT)
               DUMMY2 => DUMMY2%NEXT
               NULLIFY(DUMMY2%NEXT)
          ENDIF
          DUMMY2%I=J

          IF (I == H) THEN
               ENDREACHED=.TRUE.
               IF (ASSOCIATED(START)) THEN
                    ALLOCATE(DUMMY%NEXT)
                    DUMMY => DUMMY%NEXT
                    NULLIFY(DUMMY%NEXT)
                    DUMMY%I=I
                    DUMMY%J=TSINDEX(DUMMY%I,J)
               ELSE ! IT IS A FIRST RUN OF EXTRACT
                    ALLOCATE(START)
                    START%I=I
                    START%J=TSINDEX(START%I,J)
                    DUMMY => START
               ENDIF
               GOTO 13
          ENDIF

          E: DO K=1,MI(I)%DATA%C(1)

               DUMMY2 => START2
               IF (MI(I)%DATA%CMIN(K)==DUMMY2%I) CYCLE
               DO M=2,DEPTH
                    DUMMY2 => DUMMY2%NEXT
                    IF (MI(I)%DATA%CMIN(K)==DUMMY2%I) CYCLE E
               ENDDO

               CALL EXTRACT(H,MI(I)%DATA%CMIN(K),I)

               IF (ENDREACHED) THEN
                    ALLOCATE(DUMMY%NEXT)
                    DUMMY => DUMMY%NEXT
                    NULLIFY(DUMMY%NEXT)
                    DUMMY%I=I
                    DUMMY%J=TSINDEX(DUMMY%I,J)
                    GOTO 13
               ENDIF
          ENDDO E

     13   IF (DEPTH==1) THEN 
               NULLIFY(DUMMY2)
               DEALLOCATE(START2)
          ELSE
               DUMMY2=>START2
               DO K=2,DEPTH-1
                    DUMMY2 => DUMMY2%NEXT
               ENDDO
               DEALLOCATE(DUMMY2%NEXT)
          ENDIF
          DEPTH = DEPTH-1
     END SUBROUTINE EXTRACT

     FUNCTION TSINDEX(I,J) !RETURNS POSITION OF TS (IN ARRAY CTS OF _MINIMA I_) THAT CONNECTS I AND J 
          IMPLICIT NONE
          INTEGER :: TSINDEX
          INTEGER,INTENT(IN) :: I,J
          
          INTEGER :: K

!
! WITHOUT THIS LINE IFORT DEBUG COMPLAINS WHEN WE HIT A TSINDEX THAT IS NOT ASSIGNED
!
          TSINDEX=0
          DO K=1,MI(I)%DATA%C(1)
               IF (MI(I)%DATA%CMIN(K) == J) THEN
                    TSINDEX=K
                    RETURN
               ENDIF
          ENDDO
     END FUNCTION TSINDEX

!
! NOTE FOR STOCKT. THE PATH.<N>.XYZ FILES HAVE GOT THE CENTRE OF MASS COORDINATES PLUS A 
! UNIT VECTOR. WE NEED TO ROTATE ALL OF THESE APPROPRIATELY AS IF EACH ONE REPRESENTS
! AN ATOM, I.E. USING NEWROTGEOM, NOT NEWROTGEOMSTOCK.
!
! REMAINING PROBLEMS FOR STOCK:
! THE ENERGY IS INVARIANT TO A CHANGE IN SIGN OF ALL THE DIPOLES.
! ALIGNING PATHS PROPERLY MIGHT THEREFORE NEED AN INVERSION OF ALL
! THE DIPOLE MOMENTS. WE ARE ACTUALLY ALIGNING THE PATH.<N>.XYZ FILES,
! WHERE THE DIPOLES ARE PRESENT AS EXTRA UNIT VECTORS WHICH ARE TREATED
! AS ATOMS IN THE FOLLOWING SUBROUTINE. 
! A GOOD SOLUTION WOULD PROBABLY BE TO DEFINE AN UNAMBIGUOUS SENSE FOR THE
! DIPOLES SO THAT THIS PROBLEM DOESN'T ARISE. IT COULD ALSO CAUSE PROBLEMS
! IN RECOGNISING IDENTICAL MINIMA.
!
! MERGES PATH OUTPUT FILES TO PRODUCE FULL PATHWAY FOR THE REARRANGEMENT;
! FRAMES IN PATH ARE REVERSED AS NEEDED;
!
     SUBROUTINE MERGEXYZEOFS  
          USE KEY, ONLY: FILTH,UNRST,FILTHSTR,DEBUG,RIGIDBODY,BULKT,TWOD,STOCKT,RBAAT,PERMDIST, AMHT,SEQ,NTSITES
          USE KEYUTILS        ! FRAMES IN BITS THAT ARE GLUED TOGETHER ARE ROTATED ACCORDINGLY;
          USE COMMONS,ONLY : PARAM1,PARAM2,PARAM3
          USE AMHGLOBALS, ONLY : NMRES

          IMPLICIT NONE       ! PREREQUISITES: CHAIN OF MIN/TS CONSTRUCTED; ASSUMES PATH IS DUMPING PLUS SIDE OF THE PATH
                              ! FIRST, AND THERE ARE NO BLANK LINES AFTER LAST FRAME (!)
                              ! DOES SOMEWHAT SIMILAR WITH EOFS.TS FILES AS WELL..
          DOUBLE PRECISION RMAT(3,3), Q2(4)
          INTEGER :: I,J,K,EOF,J1,J2,INDEXTS !,FL
          DOUBLE PRECISION :: S,E,SLAST,SFIRST
          DOUBLE PRECISION,POINTER :: LASTFRAME(:)
          LOGICAL :: REVERSEORDER, EXTEST
          CHARACTER :: ITSTRING*80, EOFSSTRING*80, PSTRING*80, ESTRING*80, PSTRING2*80, RBPSTRING*80
          DOUBLE PRECISION         :: EDUMMY,TMPTS(NOPT),LGDUMMY(NOPT),RMS,QPLUS(NOPT),QMINUS(NOPT),EPLUS,EMINUS,EDUMMY2,DPRAND
          DOUBLE PRECISION         :: CMXA,CMYA,CMZA,DLENGTH,DIST,DIST2,QSAVE(NOPT),XDUMMY,X(NOPT)
          DOUBLE PRECISION         :: SITESB(3*NTSITES), SITESA(3*NTSITES), XS(3*NTSITES), CMAX,CMAY,CMAZ,CMBX,CMBY,CMBZ
          DOUBLE PRECISION         :: DISTANCE, QBEST(4), RMATBEST(3,3)
          DOUBLE PRECISION         :: CMXFIX,CMYFIX,CMZFIX,DUMMY1,GRAD(NOPT)
          LOGICAL PATHFAILT, LOCALSTOCK

!         LOCAL AMH VARIABLES
          INTEGER NRES_AMH, I_RES, GLY_COUNT, NDUMMY
          CHARACTER(LEN=5) AMHDUMMY 

          LOGICAL KNOWE, KNOWG, KNOWH
          COMMON /KNOWN/ KNOWE, KNOWG, KNOWH
          DOUBLE PRECISION BSX, BSY, BSZ
          DOUBLE PRECISION SITES(3*NTSITES)
          COMMON /BULKSHIFT/ BSX, BSY, BSZ

          TYPE LOYNO
               CHARACTER(LEN=132)       :: COMMENT
               CHARACTER(LEN=70),POINTER:: LINE(:)
               CHARACTER(LEN=5),POINTER :: SYM(:)
               DOUBLE PRECISION,POINTER          :: Q(:)
               TYPE(LOYNO),POINTER      :: NEXT
               TYPE(LOYNO),POINTER      :: PREVIOUS
          END TYPE
          TYPE(LOYNO),POINTER :: HEAD,TMP,TAIL
          NULLIFY(HEAD,TMP,TAIL,LASTFRAME)
          
          LOCALSTOCK=.FALSE.
          IF (STOCKT) LOCALSTOCK=.TRUE.
          IF (FILTH.EQ.0) THEN
             IF (RBAAT) OPEN(UNIT=55,FILE='RBPATH.XYZ',STATUS='UNKNOWN')
!            OPEN(UNIT=50,FILE='PATH.XYZ',STATUS='REPLACE',ACTION='WRITE')
             OPEN(UNIT=50,FILE='PATH.XYZ',STATUS='UNKNOWN')
             IF (UNRST) OPEN(UNIT=51,FILE='PATH.UNR.XYZ',STATUS='REPLACE',ACTION='WRITE')
!            OPEN(UNIT=41,FILE='EOFS',STATUS='REPLACE',ACTION='WRITE')
             OPEN(UNIT=41,FILE='EOFS',STATUS='UNKNOWN')
          ELSE
             IF (RBAAT) THEN
                WRITE(RBPSTRING,'(A)') 'RBPATH.XYZ.'//TRIM(ADJUSTL(FILTHSTR))
!               OPEN(UNIT=55,FILE=PSTRING,STATUS='REPLACE',ACTION='WRITE')
                OPEN(UNIT=55,FILE=PSTRING,STATUS='UNKNOWN')
             ENDIF
             WRITE(PSTRING,'(A)') 'PATH.XYZ.'//TRIM(ADJUSTL(FILTHSTR))
             IF (UNRST) WRITE(PSTRING2,'(A)') 'PATH.UNR.XYZ.'//TRIM(ADJUSTL(FILTHSTR))
             WRITE(ESTRING,'(A)') 'EOFS.'//TRIM(ADJUSTL(FILTHSTR))
!            OPEN(UNIT=50,FILE=PSTRING,STATUS='REPLACE',ACTION='WRITE')
             OPEN(UNIT=50,FILE=PSTRING,STATUS='UNKNOWN')
             IF (UNRST) OPEN(UNIT=51,FILE=PSTRING2,STATUS='REPLACE',ACTION='WRITE')
!            OPEN(UNIT=41,FILE=ESTRING,STATUS='REPLACE',ACTION='WRITE')
             OPEN(UNIT=41,FILE=ESTRING,STATUS='UNKNOWN')
          ENDIF
          
          DUMMY => START
          SLAST=0.0D0
          I=1      ! DC430 > I LOOPS OVER PATHS
          DO
             IF (.NOT.ASSOCIATED(DUMMY%NEXT)) EXIT
             REVERSEORDER=.TRUE.
             IF (TS( MI(DUMMY%I)%DATA%CTS(DUMMY%J) )%DATA%P == DUMMY%I) REVERSEORDER=.FALSE.
             CALL MKFNAMES(MI(DUMMY%I)%DATA%CTS(DUMMY%J),FILTH,FILTHSTR,ITSTRING,EOFSSTRING)
!
!  IF READSP IS TRUE THEN THE PATHWAY COORDINATES AND ENERGY FILES MAY NOT BE PRESENT.
!  WE COULD EITHER IGNORE THIS AND SKIP ON, OR ATTEMPT TO RECALCULATE THE PATHWAY.
!  IF WE RECALCULATE THE PATHWAY, THERE IS A DANGER THAT IT MIGHT NOT LINK THE EXPECTED
!  MINIMA!
!
             INQUIRE(FILE=ITSTRING,EXIST=EXTEST)
             IF (.NOT.EXTEST) THEN
                REVERSEORDER=.TRUE.
                INDEXTS=MI(DUMMY%I)%DATA%CTS(DUMMY%J)
                PRINT '(A,I6,A)',' PATHWAY FILES FOR TS ',INDEXTS,' ARE MISSING, ATTEMPTING TO RECREATE THEM'
                TMPTS(1:NOPT)=TS(INDEXTS)%DATA%X(1:NOPT)
                EDUMMY=TS(INDEXTS)%DATA%E
                ALLOCATE(TS(INDEXTS)%DATA%VECS(NOPT))
                DO J1=1,NOPT
                   TS(INDEXTS)%DATA%VECS(J1)=(DPRAND()-0.5D0)*2 ! SINCE VECS IS PRESUMABLY UNKNOWN
                   TS(INDEXTS)%DATA%EVALMIN=-1.0D0 ! THIS EIGENVALUE SHOULD BE DETERMINED
                ENDDO
                KNOWE=.FALSE.; KNOWG=.FALSE.; KNOWH=.FALSE.
                CALL PATH(TMPTS,EDUMMY,LGDUMMY,RMS,TS(INDEXTS)%DATA%EVALMIN,TS(INDEXTS)%DATA%VECS,.TRUE., &
                    & QPLUS,QMINUS,DEBUG,EDUMMY2,EPLUS,EMINUS, &
                    & TS(INDEXTS)%DATA%SLENGTH,TS(INDEXTS)%DATA%DISP,TS(INDEXTS)%DATA%GAMMA,TS(INDEXTS)%DATA%NTILDE, &
                    & FRQSTS,FRQSPLUS, FRQSMINUS,ITSTRING,EOFSSTRING,PATHFAILT)
                IF (.NOT.PATHFAILT) THEN
                   PATHFAILT=.TRUE.
                   IF (ABS(EPLUS-MI(TS(INDEXTS)%DATA%P)%DATA%E).LE.EDIFFTOL) THEN
                      PRINT '(A)',' ENERGY FROM PLUS SIDE OF PATH MATCHES EXPECTED PLUS MINIMUM'

                      CALL MINPERMDIST(MI(TS(INDEXTS)%DATA%P)%DATA%X,QPLUS,NATOMS,DEBUG, &
    &                         PARAM1,PARAM2,PARAM3,BULKT,TWOD,DIST,DIST2,RIGIDBODY,RMAT)
                      IF (D.LE.GEOMDIFFTOL) THEN
                         PRINT '(A)',' GEOMETRY FROM PLUS SIDE OF PATH MATCHES EXPECTED PLUS MINIMUM'
                         IF (ABS(EMINUS-MI(TS(INDEXTS)%DATA%M)%DATA%E).LE.EDIFFTOL) THEN
                            PRINT '(A)',' ENERGY FROM MINUS SIDE OF PATH MATCHES EXPECTED MINUS MINIMUM'
                            CALL MINPERMDIST(MI(TS(INDEXTS)%DATA%M)%DATA%X,QMINUS,NATOMS,DEBUG, &
    &                                           PARAM1,PARAM2,PARAM3,BULKT,TWOD,DIST,DIST2,RIGIDBODY,RMAT)
                            IF (D.LE.GEOMDIFFTOL) THEN
                               PRINT '(A)',' GEOMETRY FROM MINUS SIDE OF PATH MATCHES EXPECTED MINUS MINIMUM'
                               IF (TS(INDEXTS)%DATA%P == DUMMY%I) REVERSEORDER=.FALSE.
                               PATHFAILT=.FALSE.
                            ELSE
                               PRINT '(A)',' GEOMETRY FROM MINUS SIDE OF PATH DOES NOT MATCH EXPECTED MINUS MINIMUM'
                            ENDIF
                         ELSE
                            PRINT '(A)',' ENERGY FROM MINUS SIDE OF PATH DOES NOT MATCH EXPECTED MINUS MINIMUM'
                         ENDIF
                      ELSE
                         PRINT '(A)',' GEOMETRY FROM PLUS SIDE OF PATH DOES NOT MATCH THAT EXPECTED FOR THE PLUS MINIMUM'
                      ENDIF
                   ENDIF
                   IF (PATHFAILT) THEN
                      IF (ABS(EPLUS-MI(TS(INDEXTS)%DATA%M)%DATA%E).LE.EDIFFTOL) THEN
                         PRINT '(A)',' ENERGY FROM PLUS SIDE OF PATH MATCHES EXPECTED MINUS MINIMUM'
                         CALL MINPERMDIST(MI(TS(INDEXTS)%DATA%M)%DATA%X,QPLUS,NATOMS,DEBUG, &
    &                                        PARAM1,PARAM2,PARAM3,BULKT,TWOD,DIST,DIST2,RIGIDBODY,RMAT)
                         IF (D.LE.GEOMDIFFTOL) THEN
                            PRINT '(A)',' GEOMETRY FROM PLUS SIDE OF PATH MATCHES EXPECTED MINUS MINIMUM'
                            IF (ABS(EMINUS-MI(TS(INDEXTS)%DATA%P)%DATA%E).LE.EDIFFTOL) THEN
                               PRINT '(A)',' ENERGY FROM MINUS SIDE OF PATH MATCHES EXPECTED PLUS MINIMUM'
                               CALL MINPERMDIST(MI(TS(INDEXTS)%DATA%P)%DATA%X,QMINUS,NATOMS,DEBUG, &
    &                                              PARAM1,PARAM2,PARAM3,BULKT,TWOD,DIST,DIST2,RIGIDBODY,RMAT)
                               IF (D.LE.GEOMDIFFTOL) THEN
                                  PRINT '(A)',' GEOMETRY FROM MINUS SIDE OF PATH MATCHES EXPECTED PLUS MINIMUM'
                                  PATHFAILT=.FALSE.
                                  IF (TS(INDEXTS)%DATA%M == DUMMY%I) REVERSEORDER=.FALSE.
                               ELSE
                                  PRINT '(A)',' GEOMETRY FROM MINUS SIDE OF PATH DOES NOT MATCH EXPECTED PLUS MINIMUM'
                               ENDIF
                            ELSE
                               PRINT '(A)',' ENERGY FROM MINUS SIDE OF PATH DOES NOT MATCH EXPECTED PLUS MINIMUM'
                            ENDIF
                         ELSE
                            PRINT '(A)',' GEOMETRY FROM PLUS SIDE OF PATH DOES NOT MATCH THAT EXPECTED FOR THE MINUS MINIMUM'
                         ENDIF
                      ENDIF
                   ENDIF
                ENDIF
                IF (PATHFAILT) PRINT '(A)',' ERROR - ATTEMPT TO FILL IN PATHWAY IS NOT CONSISTENT WITH SAVED DATA'
             ENDIF
               
             !XYZ BIT

             !READING PART OF THE PATH THAT CORRESPONDS TO TS DUMMY%J
             OPEN(UNIT=40,FILE=ITSTRING,STATUS='OLD',ACTION="READ")
             K=1   ! DC430 > K LOOPS OVER FRAMES IN A GIVEN PATH CORRESPONDING TO A SINGLE TS
             DO
                  READ(40,*,IOSTAT=EOF) !NEVIDOMO SKIL'KY FRAIMIV
                  IF (EOF==0) THEN
                       BACKSPACE 40
                  ELSE
                       EXIT
                  ENDIF
                  IF (K<0) PRINT*,'K=',K ! STUPID FIX FOR STUPID SUN COMPILER BUG
                                         ! NEED TO ACCESS K TO PREVENT SEGV       DJW
                    
                  ! STVORYTY NOVYI ELEMENT
                  ALLOCATE(TAIL); NULLIFY(TAIL%NEXT,TAIL%PREVIOUS,TAIL%Q,TAIL%SYM,TAIL%LINE)
                    
                  ! I ZAPOVNYTY YOGO DANYMY
                  ALLOCATE(TAIL%Q(3*NATOMS),TAIL%SYM(NATOMS))
                  READ(40,'(A)')
                  READ(40,'(A)') TAIL%COMMENT
                  IF (STOCKT) THEN
                     DO J=1,(NATOMS/2)
                        READ(40,'(A5,1X,3F20.10,TR13,3F20.10)') TAIL%SYM(J), &
  &                     TAIL%Q(3*(J-1)+1),TAIL%Q(3*(J-1)+2),TAIL%Q(3*(J-1)+3), &
  &                     TAIL%Q(3*((NATOMS/2)+J-1)+1),TAIL%Q(3*((NATOMS/2)+J-1)+2),TAIL%Q(3*((NATOMS/2)+J-1)+3)
                        DLENGTH=TAIL%Q(3*((NATOMS/2)+J-1)+1)**2 &
  &                            +TAIL%Q(3*((NATOMS/2)+J-1)+2)**2 & 
  &                            +TAIL%Q(3*((NATOMS/2)+J-1)+3)**2 
                        IF (ABS(DLENGTH-1.0D0).GT.1.0D-5) THEN
                           PRINT '(A,I8,G20.10)','1 MERGE> J,LENGTH=',J,DLENGTH
                           PRINT '(A,3G20.10)','UNIT VECTOR: ',TAIL%Q(3*((NATOMS/2)+J-1)+1:3*((NATOMS/2)+J-1)+3)
                           STOP
                        ENDIF
                     ENDDO
                  ELSE IF(RBAAT) THEN
                     DO J=1,(NATOMS/2)
                        READ(40,'(A5,1X,6F20.10)') TAIL%SYM(J), &
  &                     TAIL%Q(3*(J-1)+1),TAIL%Q(3*(J-1)+2),TAIL%Q(3*(J-1)+3), &
  &                     TAIL%Q(3*((NATOMS/2)+J-1)+1),TAIL%Q(3*((NATOMS/2)+J-1)+2),TAIL%Q(3*((NATOMS/2)+J-1)+3)
                     ENDDO
                  ELSEIF (AMHT) THEN
                     GLY_COUNT = 0
                     NDUMMY=0
                     DO J=1,NMRES
                        IF (SEQ(J).EQ.8) THEN
                           READ(40,'(A5,1X,3F20.10)')AMHDUMMY,XDUMMY,XDUMMY,XDUMMY
                           NDUMMY=NDUMMY+1
                           READ(40,'(A5,1X,3F20.10)') TAIL%SYM(NDUMMY),TAIL%Q(3*(NDUMMY-1)+1),TAIL%Q(3*(NDUMMY-1)+2),  &
     &                                                TAIL%Q(3*(NDUMMY-1)+3)
                           NDUMMY=NDUMMY+1
                           READ(40,'(A5,1X,3F20.10)') TAIL%SYM(NDUMMY),TAIL%Q(3*(NDUMMY-1)+1),TAIL%Q(3*(NDUMMY-1)+2),  &
     &                                                TAIL%Q(3*(NDUMMY-1)+3)
                        ELSE
                           NDUMMY=NDUMMY+1
                           READ(40,'(A5,1X,3F20.10)') TAIL%SYM(NDUMMY),TAIL%Q(3*(NDUMMY-1)+1),TAIL%Q(3*(NDUMMY-1)+2),  &
     &                                                TAIL%Q(3*(NDUMMY-1)+3)
                           NDUMMY=NDUMMY+1
                           READ(40,'(A5,1X,3F20.10)') TAIL%SYM(NDUMMY),TAIL%Q(3*(NDUMMY-1)+1),TAIL%Q(3*(NDUMMY-1)+2),  &
     &                                                TAIL%Q(3*(NDUMMY-1)+3)
                           NDUMMY=NDUMMY+1
                           READ(40,'(A5,1X,3F20.10)') TAIL%SYM(NDUMMY),TAIL%Q(3*(NDUMMY-1)+1),TAIL%Q(3*(NDUMMY-1)+2),  &
     &                                                TAIL%Q(3*(NDUMMY-1)+3)
                        ENDIF
                     ENDDO
                  ELSE
                     DO J=1,NATOMS
                        READ(40,'(A5,1X,3F20.10)') TAIL%SYM(J),TAIL%Q(3*(J-1)+1),TAIL%Q(3*(J-1)+2),TAIL%Q(3*(J-1)+3)
                     ENDDO
                  ENDIF
                  
                  IF (ASSOCIATED(TMP)) THEN ! ODYN ABO DEKIL'KA LANOK VGE ISNYE; TEPER TMP MISTYT' POPEREDNU LANKU
                       TMP%NEXT => TAIL
                       TAIL%PREVIOUS => TMP ! ZVOROTNIY ZVYAZOK
                  ENDIF

                  TMP=>TAIL 
                  IF (K==1) HEAD => TAIL
                  K=K+1
             ENDDO
             CLOSE(40)
               
!  WRITING FRAMES TO PATH.XYZ FILE IN CORRECT ORDER + ROTATED AS NEEDED
!  FOR STOCKMAYER THE COORDINATES READ IN FROM PATH.*.XYZ CONTAIN THE
!  VECTOR COMPONENTS, RATHER THAN THAN SPHERICAL POLAR ANGLES. HENCE
!  WE WANT TO CALL NEWROTGEOM TO ROTATE ALL OF THE SITES.
!  WE ALSO NEED STOCKT TURNED OFF IN MINPERMDIST, WHICH IS WHY WE
!  USE A LOCALSTOCK IN THIS SUBROUTINE.

             K=1
             IF (REVERSEORDER) THEN
                  TMP => TAIL; NULLIFY(TAIL)
                  CMXA=0.0D0; CMYA=0.0D0; CMZA=0.0D0
                  IF (RBAAT) THEN
                     CMAX = 0.0D0; CMAY = 0.0D0; CMAZ = 0.0D0
                     DO J1 = 1, (NATOMS/2)
                        CMAX = CMAX + TMP%Q(3*(J1-1)+1)
                        CMAY = CMAY + TMP%Q(3*(J1-1)+2)
                        CMAZ = CMAZ + TMP%Q(3*(J1-1)+3)
                     ENDDO
                     CMAX = 2*CMAX/NATOMS; CMAY = 2*CMAY/NATOMS; CMAZ = 2*CMAZ/NATOMS
                  ELSE
                     DO J1=1,NATOMS
                        CMXA=CMXA+TMP%Q(3*(J1-1)+1)
                        CMYA=CMYA+TMP%Q(3*(J1-1)+2)
                        CMZA=CMZA+TMP%Q(3*(J1-1)+3)
                     ENDDO
                     CMXA=CMXA/NATOMS; CMYA=CMYA/NATOMS; CMZA=CMZA/NATOMS
                  ENDIF
                  IF (I.EQ.1) THEN
                     CMXFIX=CMXA; CMYFIX=CMYA; CMZFIX=CMZA;
                  ENDIF
                  IF (LOCALSTOCK) STOCKT=.FALSE.
                  IF ((I>1).AND.(.NOT.UNRST)) THEN
!
!  MAY NEED TO INVERT DIPOLES FOR STOCK.
!
                     IF (LOCALSTOCK) QSAVE(1:NOPT)=TMP%Q(1:NOPT)
                     IF (PERMDIST) THEN
                        IF (RBAAT) THEN
                           CALL RBMINPERMDIST(LASTFRAME,TMP%Q,DIST,DIST2,Q2,RMATBEST,DEBUG,PARAM1,PARAM2,PARAM3,BULKT,SITESB,SITESA)
                        ELSE 
                           CALL MINPERMDIST(LASTFRAME,TMP%Q,NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,DIST,DIST2,RIGIDBODY,RMAT)
                        ENDIF
                     ELSE
                        IF (RBAAT) THEN
                           CALL RBMINDIST2(LASTFRAME,TMP%Q,NATOMS,DIST,Q2,DEBUG)
                           CALL RBROT(TMP%Q,X,Q2,NATOMS)
                           DO J1 = 1, NATOMS/2
                              J2 = 3*J1
                              TMP%Q(J2-2) = X(J2-2) + CMXA
                              TMP%Q(J2-1) = X(J2-1) + CMYA
                              TMP%Q(J2)   = X(J2)   + CMZA
                           ENDDO
                        ELSE
                           CALL NEWMINDIST(LASTFRAME,TMP%Q,NATOMS,D,BULKT,TWOD,TMP%SYM(1),.TRUE.,RIGIDBODY,DEBUG,RMAT)
                           CALL NEWROTGEOM(NATOMS,TMP%Q,RMAT,CMXA,CMYA,CMZA)
                           DO J1=1,NATOMS
                              J2=3*J1
                              TMP%Q(J2-2)=TMP%Q(J2-2)+CMXFIX-CMXA
                              TMP%Q(J2-1)=TMP%Q(J2-1)+CMYFIX-CMYA
                              TMP%Q(J2)  =TMP%Q(J2)  +CMZFIX-CMZA
                           ENDDO
                        ENDIF
                     ENDIF
                  ENDIF
                  IF (LOCALSTOCK) STOCKT=.TRUE.
                  DO
                     IF (RBAAT) THEN
                        CMAX = 0.0D0; CMAY = 0.0D0; CMAZ = 0.0D0
                        DO J1 = 1, (NATOMS/2)
                           CMAX = CMAX + TMP%Q(3*(J1-1)+1)
                           CMAY = CMAY + TMP%Q(3*(J1-1)+2)
                           CMAZ = CMAZ + TMP%Q(3*(J1-1)+3)
                        ENDDO
                        CMAX = 2*CMAX/NATOMS; CMAY = 2*CMAY/NATOMS; CMAZ = 2*CMAZ/NATOMS
                        DO J1 = 1, NATOMS/2
                           TMP%Q(3*(J1-1)+1) =-CMAX + TMP%Q(3*(J1-1)+1)
                           TMP%Q(3*(J1-1)+2) =-CMAY + TMP%Q(3*(J1-1)+2)
                           TMP%Q(3*(J1-1)+3) =-CMAZ + TMP%Q(3*(J1-1)+3)
                        ENDDO
                        IF (I == 1) THEN
                           CALL SITEPOS(TMP%Q,XS)
                        ELSEIF (K == 1) THEN
                           XS(:) = SITESA(:)
                           DO J1 = 1, NATOMS/2
                              J2 = 3*J1
                              TMP%Q(J2-2:J2) = MATMUL(RMATBEST,TMP%Q(J2-2:J2))
                           ENDDO
                        ELSEIF (I > 1 .AND. K > 1) THEN
                           CALL SITEPOS(TMP%Q,XS)
                           DO J1 = 1, NTSITES
                              J2 = 3*J1
                              XS(J2-2:J2) = MATMUL(RMATBEST,XS(J2-2:J2))
                           ENDDO
                           DO J1 = 1, NATOMS/2
                              J2 = 3*J1
                              TMP%Q(J2-2:J2) = MATMUL(RMATBEST,TMP%Q(J2-2:J2))
                           ENDDO
                        ENDIF
                        CALL RBFRAME(TMP%COMMENT,XS,TMP%Q)
                     ELSE
                        IF ((I>1.AND.K>1) .AND. (.NOT.UNRST)) THEN
                           IF (.NOT.BULKT) THEN
                              CALL NEWROTGEOM(NATOMS,TMP%Q,RMAT,CMXA,CMYA,CMZA)
                              DO J1=1,NATOMS
                                 J2=3*J1
                                 TMP%Q(J2-2)=TMP%Q(J2-2)+CMXFIX-CMXA
                                 TMP%Q(J2-1)=TMP%Q(J2-1)+CMYFIX-CMYA
                                 TMP%Q(J2)  =TMP%Q(J2)  +CMZFIX-CMZA
                              ENDDO
                           ELSE
                              DO J2=1,NATOMS
                                 TMP%Q(3*(J2-1)+1)=TMP%Q(3*(J2-1)+1)-BSX-PARAM1*NINT((TMP%Q(3*(J2-1)+1)-BSX)/PARAM1)
                                 TMP%Q(3*(J2-1)+2)=TMP%Q(3*(J2-1)+2)-BSY-PARAM2*NINT((TMP%Q(3*(J2-1)+2)-BSY)/PARAM2)
                                 IF (.NOT.TWOD) TMP%Q(3*(J2-1)+3)= &
  &                                 TMP%Q(3*(J2-1)+3)-BSZ-PARAM3*NINT((TMP%Q(3*(J2-1)+3)-BSZ)/PARAM3)
                              ENDDO
                           ENDIF
                        ENDIF
                        CALL WRITEFRAME(TMP%COMMENT,TMP%SYM,TMP%Q)
                     ENDIF
                     IF (UNRST) CALL WRITEFRAMEUNRES(TMP%COMMENT,TMP%SYM,TMP%Q)
                     IF (.NOT.ASSOCIATED(TMP%PREVIOUS)) THEN
                          NULLIFY(HEAD)
                          IF (I>1) DEALLOCATE(LASTFRAME)
                          LASTFRAME => TMP%Q
                          DEALLOCATE(TMP%SYM)
                          DEALLOCATE(TMP)
                          EXIT
                          
                     ELSE
                          DEALLOCATE(TMP%Q,TMP%SYM)
                          TMP => TMP%PREVIOUS
                          NULLIFY(TMP%NEXT%PREVIOUS)
                          DEALLOCATE(TMP%NEXT)
                     ENDIF
                     K=K+1
                  ENDDO
             ELSE
                  TMP => HEAD; NULLIFY(HEAD)
                  CMXA=0.0D0; CMYA=0.0D0; CMZA=0.0D0
                  IF (RBAAT) THEN
                     CMAX = 0.0D0; CMAY = 0.0D0; CMAZ = 0.0D0
                     DO J1 = 1, (NATOMS/2)
                        CMAX = CMAX + TMP%Q(3*(J1-1)+1)
                        CMAY = CMAY + TMP%Q(3*(J1-1)+2)
                        CMAZ = CMAZ + TMP%Q(3*(J1-1)+3)
                     ENDDO
                     CMAX = 2*CMAX/NATOMS; CMAY = 2*CMAY/NATOMS; CMAZ = 2*CMAZ/NATOMS
                  ELSE
                     DO J1=1,NATOMS
                        CMXA=CMXA+TMP%Q(3*(J1-1)+1)
                        CMYA=CMYA+TMP%Q(3*(J1-1)+2)
                        CMZA=CMZA+TMP%Q(3*(J1-1)+3)
                     ENDDO
                     CMXA=CMXA/NATOMS; CMYA=CMYA/NATOMS; CMZA=CMZA/NATOMS
                  ENDIF
                  IF (I.EQ.1) THEN
                     CMXFIX=CMXA; CMYFIX=CMYA; CMZFIX=CMZA;
                  ENDIF
                  IF ((I>1).AND.(.NOT.UNRST)) THEN
                     IF (LOCALSTOCK) STOCKT=.FALSE.
                     IF (PERMDIST) THEN
                        IF (RBAAT) THEN
                           CALL RBMINPERMDIST(LASTFRAME,TMP%Q,DIST,DIST2,Q2,RMATBEST,DEBUG,PARAM1,PARAM2,PARAM3,BULKT,SITESB,SITESA)
                        ELSE
                           CALL MINPERMDIST(LASTFRAME,TMP%Q,NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,DIST,DIST2,RIGIDBODY,RMAT)
                        ENDIF
                     ELSE
                        IF (RBAAT) THEN
                           CALL RBMINDIST2(LASTFRAME,TMP%Q,NATOMS,DIST,Q2,DEBUG)
                           CALL RBROT(TMP%Q,X,Q2,NATOMS)
                           DO J1 = 1, NATOMS/2
                              J2 = 3*J1
                              TMP%Q(J2-2) = X(J2-2) + CMXA
                              TMP%Q(J2-1) = X(J2-1) + CMYA
                              TMP%Q(J2)   = X(J2)   + CMZA
                           ENDDO
                        ELSE
                           CALL NEWMINDIST(LASTFRAME,TMP%Q,NATOMS,D,BULKT,TWOD,TMP%SYM(1),.TRUE.,RIGIDBODY,DEBUG,RMAT)
                           CALL NEWROTGEOM(NATOMS,TMP%Q,RMAT,CMXA,CMYA,CMZA)
                           DO J1=1,NATOMS
                              J2=3*J1
                              TMP%Q(J2-2)=TMP%Q(J2-2)+CMXFIX-CMXA
                              TMP%Q(J2-1)=TMP%Q(J2-1)+CMYFIX-CMYA
                              TMP%Q(J2)  =TMP%Q(J2)  +CMZFIX-CMZA
                           ENDDO
                        ENDIF
                     ENDIF
                     IF (LOCALSTOCK) STOCKT=.TRUE.
                  ENDIF
                  DO
                     IF (RBAAT) THEN
                        CMAX = 0.0D0; CMAY = 0.0D0; CMAZ = 0.0D0
                        DO J1 = 1, (NATOMS/2)
                           CMAX = CMAX + TMP%Q(3*(J1-1)+1)
                           CMAY = CMAY + TMP%Q(3*(J1-1)+2)
                           CMAZ = CMAZ + TMP%Q(3*(J1-1)+3)
                        ENDDO
                        CMAX = 2*CMAX/NATOMS; CMAY = 2*CMAY/NATOMS; CMAZ = 2*CMAZ/NATOMS
                        DO J1 = 1, NATOMS/2
                           TMP%Q(3*(J1-1)+1) =-CMAX + TMP%Q(3*(J1-1)+1)
                           TMP%Q(3*(J1-1)+2) =-CMAY + TMP%Q(3*(J1-1)+2)
                           TMP%Q(3*(J1-1)+3) =-CMAZ + TMP%Q(3*(J1-1)+3)
                        ENDDO
                        IF (I == 1) THEN
                           CALL SITEPOS(TMP%Q,XS)
                        ELSEIF (K == 1) THEN
                           XS(:) = SITESA(:)
                           DO J1 = 1, NATOMS/2
                              J2 = 3*J1
                              TMP%Q(J2-2:J2) = MATMUL(RMATBEST,TMP%Q(J2-2:J2))
                           ENDDO
                        ELSEIF (I > 1 .AND. K > 1) THEN
                           CALL SITEPOS(TMP%Q,XS)
                           DO J1 = 1, NTSITES
                              J2 = 3*J1
                              XS(J2-2:J2) = MATMUL(RMATBEST,XS(J2-2:J2))
                           ENDDO
                           DO J1 = 1, NATOMS/2
                              J2 = 3*J1
                              TMP%Q(J2-2:J2) = MATMUL(RMATBEST,TMP%Q(J2-2:J2))
                           ENDDO
                        ENDIF
                        CALL RBFRAME(TMP%COMMENT,XS,TMP%Q)
                     ELSE
                        IF ((I>1.AND.K>1).AND.(.NOT.UNRST)) THEN
                           IF (.NOT.BULKT) THEN
                              CALL NEWROTGEOM(NATOMS,TMP%Q,RMAT,CMXA,CMYA,CMZA) 
                              DO J1=1,NATOMS
                                 TMP%Q(3*(J1-1)+1)=TMP%Q(3*(J1-1)+1)-CMXA+CMXFIX
                                 TMP%Q(3*(J1-1)+2)=TMP%Q(3*(J1-1)+2)-CMYA+CMYFIX
                                 TMP%Q(3*(J1-1)+3)=TMP%Q(3*(J1-1)+3)-CMZA+CMZFIX
                              ENDDO
                           ELSE
                              DO J2=1,NATOMS
                                 TMP%Q(3*(J2-1)+1)=TMP%Q(3*(J2-1)+1)-BSX-PARAM1*NINT((TMP%Q(3*(J2-1)+1)-BSX)/PARAM1)
                                 TMP%Q(3*(J2-1)+2)=TMP%Q(3*(J2-1)+2)-BSY-PARAM2*NINT((TMP%Q(3*(J2-1)+2)-BSY)/PARAM2)
                                 IF (.NOT.TWOD) TMP%Q(3*(J2-1)+3)= &
     &                                               TMP%Q(3*(J2-1)+3)-BSZ-PARAM3*NINT((TMP%Q(3*(J2-1)+3)-BSZ)/PARAM3)
                              ENDDO
                           ENDIF
                        ENDIF
                        CALL WRITEFRAME(TMP%COMMENT,TMP%SYM,TMP%Q)
                     ENDIF
                     IF (UNRST) CALL WRITEFRAMEUNRES(TMP%COMMENT,TMP%SYM,TMP%Q)
                     IF (.NOT.ASSOCIATED(TMP%NEXT)) THEN
                        NULLIFY(TAIL)
                        IF (I>1) DEALLOCATE(LASTFRAME)
                           LASTFRAME => TMP%Q
                           DEALLOCATE(TMP%SYM)
                           DEALLOCATE(TMP)
                        EXIT
                     ELSE
                        DEALLOCATE(TMP%Q,TMP%SYM)
                        TMP => TMP%NEXT
                        NULLIFY(TMP%PREVIOUS%NEXT)
                        DEALLOCATE(TMP%PREVIOUS)
                     ENDIF
                     K=K+1
                  ENDDO
             ENDIF

             !ENERGY BIT
             IF (REVERSEORDER) THEN
                  !OPEN(UNIT=40,FILE=EOFSSTRING,STATUS='OLD',ACTION='READ',POSITION='APPEND')
                  ! PGI COMPLAINS ABOUT READ AND APPEND TOGETHER ... HM...
                  OPEN(UNIT=40,FILE=EOFSSTRING,STATUS='OLD',POSITION='APPEND')
                  BACKSPACE 40
                  READ(40,'(F20.10)') SFIRST
                  REWIND 40
                  SFIRST=SFIRST+SLAST
                  DO
                       READ(40,*,IOSTAT=EOF) !TRY TO READ
                       IF (EOF==0) THEN
                            BACKSPACE 40
                       ELSE
                            EXIT
                       ENDIF

                       ALLOCATE(HEAD)
                       NULLIFY(HEAD%NEXT)
                       IF (ASSOCIATED(TMP)) HEAD%NEXT => TMP
                       ALLOCATE(HEAD%LINE(1))
                       READ(40,'(A)') HEAD%LINE(1)
                       TMP=>HEAD
                  ENDDO
                  DO
                       READ(HEAD%LINE(1),'(2F20.10)') S,E
                       DEALLOCATE(HEAD%LINE)
                       WRITE(UNIT=41,FMT='(2F20.10)') SFIRST - S,E
                       SLAST = SFIRST - S
                       IF (.NOT.ASSOCIATED(HEAD%NEXT)) EXIT
                       TMP => HEAD%NEXT
                       NULLIFY(HEAD%NEXT)
                       DEALLOCATE(HEAD)
                       HEAD => TMP
                  ENDDO
                  NULLIFY(HEAD)
                  DEALLOCATE(TMP)
                  NULLIFY(TMP)
             ELSE
                  OPEN(UNIT=40,FILE=EOFSSTRING,STATUS='OLD',ACTION='READ')
                  J=1
                  DO
                     READ(UNIT=40,FMT='(2F20.10)',IOSTAT=EOF) S,E
                     IF (.NOT.EOF==0) EXIT
                       IF (J==1) SFIRST=ABS(S) + SLAST
                       WRITE(UNIT=41,FMT='(2F20.10)') SFIRST + S,E
                       SLAST=SFIRST + S
                       J=0
                  ENDDO
             ENDIF

             CLOSE(40)
               
             DUMMY=>DUMMY%NEXT
             I=I+1
          ENDDO
          DEALLOCATE(LASTFRAME)
          CLOSE(50)
          IF (RBAAT) CLOSE(55)
          CLOSE(41)
     END SUBROUTINE MERGEXYZEOFS

     SUBROUTINE WRITEFRAME(C,S,Q)
          USE KEY,ONLY : STOCKT, RBAAT, NTSITES
          IMPLICIT NONE

          CHARACTER(LEN=132),INTENT(IN)         :: C
          CHARACTER(LEN=5),POINTER,DIMENSION(:) :: S
          DOUBLE PRECISION,POINTER,DIMENSION(:) :: Q
          DOUBLE PRECISION SITES(3*NTSITES)

          INTEGER :: J

          IF (STOCKT) THEN
             WRITE(50,'(I6)') (NATOMS/2)
             WRITE(50,'(1X,A)') TRIM(ADJUSTL(C))
             DO J=1,(NATOMS/2)
!               WRITE(50,'(A5,1X,3F20.10,A13,3F20.10)') S(J),Q(3*(J-1)+1),Q(3*(J-1)+2),Q(3*(J-1)+3), &
! &                             ' ATOM_VECTOR ', &
! &                             Q(3*((NATOMS/2)+J-1)+1),Q(3*((NATOMS/2)+J-1)+2),Q(3*((NATOMS/2)+J-1)+3)
                WRITE(50,'(A5,1X,6F20.10)') S(J),Q(3*(J-1)+1),Q(3*(J-1)+2),Q(3*(J-1)+3), &
  &                             Q(3*((NATOMS/2)+J-1)+1),Q(3*((NATOMS/2)+J-1)+2),Q(3*((NATOMS/2)+J-1)+3)
             ENDDO
          ELSE IF (RBAAT) THEN
             WRITE(50,'(I6)') (NATOMS/2)
             WRITE(50,'(1X,A)') TRIM(ADJUSTL(C))
             WRITE(55,'(I6)') NTSITES
             WRITE(55,'(1X,A)') TRIM(ADJUSTL(C))
             DO J=1,(NATOMS/2)
                WRITE(50,'(A5,1X,6F20.10)') S(J),Q(3*(J-1)+1),Q(3*(J-1)+2),Q(3*(J-1)+3), &
  &                             Q(3*((NATOMS/2)+J-1)+1),Q(3*((NATOMS/2)+J-1)+2),Q(3*((NATOMS/2)+J-1)+3)
             ENDDO
             CALL SITEPOS(Q,SITES)
             DO J=1,NTSITES
                WRITE(55,'(A5,1X,3F20.10)') 'LA ',SITES(3*(J-1)+1),SITES(3*(J-1)+2),SITES(3*(J-1)+3)
             ENDDO
          ELSE
             WRITE(50,'(I6)') NATOMS
             WRITE(50,'(1X,A)') TRIM(ADJUSTL(C))
             DO J=1,NATOMS
                WRITE(50,'(A5,1X,3F20.10)') S(J),Q(3*(J-1)+1),Q(3*(J-1)+2),Q(3*(J-1)+3)
             ENDDO
          ENDIF
     END SUBROUTINE WRITEFRAME

     SUBROUTINE WRITEFRAMEUNRES(C,S,Q) ! JMC WRITE COORDS PLUS DUMMY PEPTIDE GROUPS FOR UNRES VISUALISATION PURPOSES
          IMPLICIT NONE

          CHARACTER(LEN=132),INTENT(IN)         :: C
          CHARACTER(LEN=5),POINTER,DIMENSION(:) :: S
          DOUBLE PRECISION,POINTER,DIMENSION(:)          :: Q
          DOUBLE PRECISION                               :: PEPCOORDS(3*NATOMS)

          INTEGER :: J,K1,K2

          WRITE(51,'(I6)') 2*NATOMS-2 ! JMC FOR CARBONS PLUS DUMMY PEPTIDE ATOMS, WANT 2*NATOMS-2...
          WRITE(51,'(1X,A)') TRIM(ADJUSTL(C))
          DO J=1,NATOMS
               WRITE(51,'(A5,1X,3F20.10)') S(J),Q(3*(J-1)+1),Q(3*(J-1)+2),Q(3*(J-1)+3)
          ENDDO
          
          DO K1=1,(NATOMS/2)-1
             DO K2=1,3
                PEPCOORDS(6*(K1-1)+K2)=(2.0D0*Q(6*(K1-1)+K2)+Q(6*K1+K2))/3.0D0
                PEPCOORDS(6*(K1-1)+K2+3)=(Q(6*(K1-1)+K2)+2.0D0*Q(6*K1+K2))/3.0D0
             ENDDO
          ENDDO
          WRITE(51,'(A2,4X,3F20.10)') ('O ',PEPCOORDS(6*(K1-1)+1),PEPCOORDS(6*(K1-1)+2),PEPCOORDS(6*(K1-1)+3) &
                                       & ,K1=1,(NATOMS/2)-1)
          WRITE(51,'(A2,4X,3F20.10)') ('N ',PEPCOORDS(6*(K1-1)+4),PEPCOORDS(6*(K1-1)+5),PEPCOORDS(6*(K1-1)+6) &
                                       & ,K1=1,(NATOMS/2)-1)

     END SUBROUTINE WRITEFRAMEUNRES

     SUBROUTINE REALLOCATETSRACK  ! RESIZES TS RACK: NEWSIZE = OLDSIZE*INCR
          IMPLICIT NONE
          
          INTEGER, PARAMETER :: INCR=2
          INTEGER :: I
          
          TEMPTSRACK => TS
          TSRACKSIZE=TSRACKSIZE*INCR
          WRITE(CHR,'(I5)') TSRACKSIZE
          PRINT '(1X,A)', 'TS RACK REDIMENTIONED TO HOLD '//TRIM(ADJUSTL(CHR))//' TS.'
          ALLOCATE(TS(TSRACKSIZE))
          DO I=1,TSRACKSIZE/INCR
               TS(I)=TEMPTSRACK(I)
          ENDDO
          DEALLOCATE(TEMPTSRACK)
     END SUBROUTINE REALLOCATETSRACK

     SUBROUTINE REALLOCATEMINRACK  ! RESIZES MINIMA RACK: NEWSIZE = OLDSIZE*INCR
          IMPLICIT NONE

          INTEGER, PARAMETER :: INCR=2
          INTEGER :: I

          TEMPMINRACK => MI
          MINRACKSIZE=MINRACKSIZE*INCR
          WRITE(CHR,'(I5)') MINRACKSIZE
          PRINT '(1X,A)', 'MINIMA RACK REDIMENTIONED TO HOLD '//TRIM(ADJUSTL(CHR))//' MINIMA.'
          ALLOCATE(MI(MINRACKSIZE))
          DO I=1,MINRACKSIZE/INCR
               MI(I)=TEMPMINRACK(I)
          ENDDO
          DEALLOCATE(TEMPMINRACK)
     END SUBROUTINE REALLOCATEMINRACK

     SUBROUTINE DUMPDB(FINISHED,FINALPATHTS)
          IMPLICIT NONE
     
          LOGICAL,INTENT(IN) :: FINISHED
          INTEGER,POINTER :: FINALPATHTS(:)

          INTEGER :: RECLEN, I,J, NSP,NMINSAVED,NTSSAVED,REC1,REC2
          INTEGER,POINTER :: MINRECORDS(:)
          LOGICAL,POINTER :: MINSAVED(:)
          CHARACTER(LEN=256) :: MYSTR
          
          ALLOCATE(MINRECORDS(NMIN),MINSAVED(NMIN))
          MINSAVED=.FALSE.
          MINRECORDS=0
          INQUIRE(IOLENGTH=RECLEN) (MI(DUMMY%I)%DATA%X(I),I=1,3*NATOMS)
          
          IF (FINISHED) THEN
               DUMMY=>START
               NSP=0
               DO
                    NSP=NSP+1
                    WRITE(MYSTR,*) NSP
                    MYSTR='POINTS'//TRIM(ADJUSTL(MYSTR))//'.OUT'
                    OPEN(UNIT=38,FILE=TRIM(ADJUSTL(MYSTR)),STATUS='UNKNOWN',FORM='UNFORMATTED',ACCESS='DIRECT',RECL=RECLEN)
                    WRITE(38,REC=1) (MI(DUMMY%I)%DATA%X(I),I=1,3*NATOMS)
                    CLOSE(38)
                    IF (ASSOCIATED(DUMMY%NEXT)) THEN
                         NSP=NSP+1
                         WRITE(MYSTR,*) NSP
                         MYSTR='POINTS'//TRIM(ADJUSTL(MYSTR))//'.OUT'
                         OPEN(UNIT=38,FILE=TRIM(ADJUSTL(MYSTR)),STATUS='UNKNOWN',FORM='UNFORMATTED',ACCESS='DIRECT',RECL=RECLEN)
                         WRITE(38,REC=1) (TS(MI(DUMMY%I)%DATA%CTS(DUMMY%J))%DATA%X(I),I=1,NOPT)
                         CLOSE(38)
                         DUMMY=>DUMMY%NEXT
                    ELSE
                         EXIT
                    ENDIF
               ENDDO
          ENDIF

          IF (ASSOCIATED(FINALPATHTS).AND.NTS==SIZE(FINALPATHTS)) THEN
               RETURN
          ELSE
               OPEN(UNIT=38,FILE="POINTS.TS",STATUS='UNKNOWN',FORM='UNFORMATTED',ACCESS='DIRECT',RECL=RECLEN)
               OPEN(UNIT=39,FILE="TS.DATA",STATUS='UNKNOWN',FORM='FORMATTED')
               OPEN(UNIT=40,FILE="POINTS.MIN",STATUS='UNKNOWN',FORM='UNFORMATTED',ACCESS='DIRECT',RECL=RECLEN)
               NMINSAVED=0
               NTSSAVED=0
               MAIN: DO J=1,NTS
                    IF (TS(J)%DATA%BAD) CYCLE ! DATA%P AND DATA%M WILL BE UNDEFINED!
                    IF (ASSOCIATED(FINALPATHTS)) THEN
                         DO I=1,SIZE(FINALPATHTS)
                              IF (FINALPATHTS(I)==J) CYCLE MAIN
                         ENDDO     
                    ENDIF
                    NTSSAVED=NTSSAVED+1
                    WRITE(38,REC=NTSSAVED) ( TS(J)%DATA%X(I),I=1,3*NATOMS )
                    IF (MINSAVED(TS(J)%DATA%P)) THEN
                         REC1=MINRECORDS(TS(J)%DATA%P)
                    ELSE
                         NMINSAVED=NMINSAVED+1
                         MINSAVED(TS(J)%DATA%P)=.TRUE.
                         MINRECORDS(TS(J)%DATA%P)=NMINSAVED
                         WRITE(40,REC=NMINSAVED) ( MI(TS(J)%DATA%P)%DATA%X(I),I=1,3*NATOMS )
                         REC1=NMINSAVED
                    ENDIF
                    IF (MINSAVED(TS(J)%DATA%M)) THEN
                         REC2=MINRECORDS(TS(J)%DATA%M)
                    ELSE
                         NMINSAVED=NMINSAVED+1
                         MINSAVED(TS(J)%DATA%M)=.TRUE.
                         MINRECORDS(TS(J)%DATA%M)=NMINSAVED
                         WRITE(40,REC=NMINSAVED) ( MI(TS(J)%DATA%M)%DATA%X(I),I=1,3*NATOMS )
                         REC2=NMINSAVED
                    ENDIF
                    WRITE(39,'(2I10)') REC1,REC2
               ENDDO MAIN
               CLOSE(39)
               CLOSE(38)
               CLOSE(40)
               IF (ASSOCIATED(FINALPATHTS)) DEALLOCATE(FINALPATHTS)
          ENDIF
     END SUBROUTINE DUMPDB

     SUBROUTINE MAKEPATHINFO
     USE SYMINF
     USE MODHESS
     USE MODCHARMM
     USE PORFUNCS
     USE MODUNRES
     USE KEY, ONLY: FILTH, FILTHSTR, UNRST, TWOD, BULKT, MACHINE, DEBUG, RIGIDBODY, NOFRQS, PERMDIST, AMHT, SDT
!    USE AMH
     USE COMMONS, ONLY: ATMASS, NINTS, ZSYM, PARAM1, PARAM2, PARAM3
     IMPLICIT NONE
     DOUBLE PRECISION RMAT(3,3), DIST, DIST2

!    LOCAL AMH VARIABLES
     INTEGER :: I_RES, GLY_COUNT
     CHARACTER(LEN=5) :: TARFL

     CHARACTER(LEN=20) :: PINFOSTRING
     DOUBLE PRECISION :: DIHE,ALLANG,DISTPF,DUMMY1,GRAD(3*NATOMS),RMS,DIAG(3*NATOMS),TEMPA(9*NATOMS),DUMQ(3*NATOMS)
     DOUBLE PRECISION :: PREVIOUSTS(3*NATOMS), INERTIA(3,3)
     INTEGER :: HORDER,I,INFO,J2,K1,RECLEN,ISTAT
     LOGICAL :: BTEST,KD,NNZ,NINTB,TSFRQDONE,MINFRQDONE

!   
!  FILE PATH.INFO IS NOW OPENED IN KEYWORD.F
!  
!    IF (FILTH.EQ.0) THEN
!       WRITE(PINFOSTRING,'(A9)') 'PATH.INFO'
!    ELSE
!       WRITE(PINFOSTRING,'(A)') 'PATH.INFO.'//TRIM(ADJUSTL(FILTHSTR))
!    ENDIF
!    IF (MACHINE) THEN
!         OPEN(UNIT=88,FILE=PINFOSTRING,STATUS='UNKNOWN',FORM='UNFORMATTED')
!    ELSE
!         OPEN(UNIT=88,FILE=PINFOSTRING,STATUS='UNKNOWN')
!    ENDIF

     TSFRQDONE=.FALSE.
     MINFRQDONE=.FALSE.

     DUMMY=>START
     I=1
     DO
        IF (UNRST.AND.CALCDIHE) THEN
           CALL UNRESCALCDIHEREF(DIHE,ALLANG,MI(DUMMY%I)%DATA%X)
!               WRITE(88,'(3F20.10)') MI(DUMMY%I)%DATA%E, DIHE, ALLANG
           WRITE(88,'(3F20.10)') MI(DUMMY%I)%DATA%E, DIHE
        ELSE
           WRITE(88,'(F20.10)') MI(DUMMY%I)%DATA%E
        ENDIF
        IF ((.NOT.MINFRQDONE).OR.(TWOD).OR.CHRMMT.OR.UNRST) THEN
           IF (UNRST) THEN ! JMC UPDATE COORDS
              DO K1=1,NRES
                  C(1,K1)=MI(DUMMY%I)%DATA%X(6*(K1-1)+1)
                  C(2,K1)=MI(DUMMY%I)%DATA%X(6*(K1-1)+2)
                  C(3,K1)=MI(DUMMY%I)%DATA%X(6*(K1-1)+3)
                  C(1,K1+NRES)=MI(DUMMY%I)%DATA%X(6*(K1-1)+4)
                  C(2,K1+NRES)=MI(DUMMY%I)%DATA%X(6*(K1-1)+5)
                  C(3,K1+NRES)=MI(DUMMY%I)%DATA%X(6*(K1-1)+6)
               ENDDO
               CALL UPDATEDC
               CALL INT_FROM_CART(.TRUE.,.FALSE.)
               CALL CHAINBUILD
           ENDIF
!          CALL POTENTIAL(MI(DUMMY%I)%DATA%X,DUMMY1,GRAD,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
           IF (CHRMMT) THEN
              HORDER=1
              FPGRP='C1'
              IF (MACHINE) THEN
                 WRITE(88) HORDER,FPGRP
              ELSE
                 WRITE(88,'(I6,1X,A4)') HORDER,FPGRP
              ENDIF
              IF (.NOT.NOFRQS) THEN
                 IF (ENDNUMHESS) THEN
                    CALL MAKENUMHESS(MI(DUMMY%I)%DATA%X,NATOMS)
                 ELSE
                    CALL POTENTIAL(MI(DUMMY%I)%DATA%X,DUMMY1,GRAD,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
                 ENDIF
                 CALL MASSWT2(NATOMS,ATMASS,DUMQ,GRAD,.TRUE.)
                 CALL DSYEV('N','U',3*NATOMS,HESS,SIZE(HESS,1),DIAG,TEMPA,9*NATOMS,INFO)
                 IF (DIAG(1).LT.DIAG(3*NATOMS)) CALL EIGENSORT_VAL_ASC(DIAG,HESS,3*NATOMS,3*NATOMS)
                 IF (MACHINE) THEN
                      WRITE(88) (DIAG(J2)*4.184D26,J2=1,3*NATOMS)
                 ELSE
                      WRITE(88,'(3G20.10)') (DIAG(J2)*4.184D26,J2=1,3*NATOMS)
                 ENDIF
              ENDIF
           ELSEIF (UNRST) THEN
              HORDER=1
              FPGRP='C1'
              WRITE(88,'(I6,1X,A4)') HORDER,FPGRP
              IF (.NOT.NOFRQS) THEN
                 IF (ENDNUMHESS) THEN
                    CALL MAKENUMINTHESS(NINTS,NATOMS)
                    CALL GETSTUFF(KD,NNZ,NINTB)
                    CALL INTSECDET(MI(DUMMY%I)%DATA%X,3*NATOMS,KD,NNZ,NINTB,DIAG)
                 ELSE
                    CALL POTENTIAL(MI(DUMMY%I)%DATA%X,DUMMY1,GRAD,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
                 ENDIF
                 WRITE(88,'(3G20.10)') (DIAG(J2),J2=1,3*NATOMS)
              ENDIF
           ELSEIF (AMHT) THEN
              IF (.NOT.NOFRQS) THEN
                 IF (ENDNUMHESS) THEN
                    CALL MAKENUMHESS(MI(DUMMY%I)%DATA%X,NATOMS)
                 ELSE
                    CALL POTENTIAL(MI(DUMMY%I)%DATA%X,DUMMY1,GRAD,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
                 ENDIF
                 CALL MASSWT(NATOMS,ATMASS,DUMQ,GRAD,.TRUE.)
                 CALL DSYEV('N','U',3*NATOMS,HESS,SIZE(HESS,1),DIAG,TEMPA,9*NATOMS,INFO)
                 IF (DIAG(1).LT.DIAG(3*NATOMS)) CALL EIGENSORT_VAL_ASC(DIAG,HESS,3*NATOMS,3*NATOMS)
                 WRITE(88,'(3G20.10)') (DIAG(J2),J2=1,3*NATOMS)
              ENDIF
           ELSE
              CALL SYMMETRY(HORDER,.FALSE.,MI(DUMMY%I)%DATA%X,INERTIA)
              WRITE(88,'(I6,1X,A4)') HORDER,FPGRP
              IF (.NOT.NOFRQS) THEN
                 IF (ENDNUMHESS) THEN
                    CALL MAKENUMHESS(MI(DUMMY%I)%DATA%X,NATOMS)
                 ELSE
                    CALL POTENTIAL(MI(DUMMY%I)%DATA%X,DUMMY1,GRAD,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
                 ENDIF
                 CALL MASSWT(NATOMS,ATMASS,DUMQ,GRAD,.TRUE.)
                 CALL DSYEV('N','U',3*NATOMS,HESS,SIZE(HESS,1),DIAG,TEMPA,9*NATOMS,INFO)
                 IF (DIAG(1).LT.DIAG(3*NATOMS)) CALL EIGENSORT_VAL_ASC(DIAG,HESS,3*NATOMS,3*NATOMS)
                 IF (SDT) THEN
                    WRITE(88,'(3G20.10)') (DIAG(J2)*4.184D26,J2=1,3*NATOMS)
                 ELSE
                    WRITE(88,'(3G20.10)') (DIAG(J2),J2=1,3*NATOMS)
                 ENDIF
              ENDIF
           ENDIF
        ELSE
           CALL SYMMETRY(HORDER,.FALSE.,MI(DUMMY%I)%DATA%X,INERTIA)
           WRITE(88,'(I6,1X,A4)') HORDER,FPGRP
! JMC           WRITE(88,'(3G20.10)') (FSAVEMIN(J2,J1),J2=1,3*NATOMS)
        ENDIF
        IF (I.GT.1) THEN
           CALL MINPERMDIST(PREVIOUSTS,MI(DUMMY%I)%DATA%X,NATOMS,DEBUG, &
    &                       PARAM1,PARAM2,PARAM3,BULKT,TWOD,DIST,DIST2,RIGIDBODY,RMAT)
           DISTPF=DIST
        ENDIF
        IF (MACHINE) THEN
             WRITE(88) (MI(DUMMY%I)%DATA%X(J2),J2=1,3*NATOMS)
        ELSE
             WRITE(88,'(3F25.15)') (MI(DUMMY%I)%DATA%X(J2),J2=1,3*NATOMS)
        ENDIF

        IF (.NOT.ASSOCIATED(DUMMY%NEXT)) EXIT

        IF (MACHINE) THEN
           WRITE(88) TS(MI(DUMMY%I)%DATA%CTS(DUMMY%J))%DATA%E
        ELSE
           WRITE(88,'(F20.10)') TS(MI(DUMMY%I)%DATA%CTS(DUMMY%J))%DATA%E
        ENDIF
        IF ((.NOT.TSFRQDONE).OR.(TWOD).OR.CHRMMT.OR.UNRST) THEN
           IF (UNRST) THEN
             DO K1=1,NRES ! JMC UPDATE COORDS
                 C(1,K1)=TS(MI(DUMMY%I)%DATA%CTS(DUMMY%J))%DATA%X(6*(K1-1)+1)
                 C(2,K1)=TS(MI(DUMMY%I)%DATA%CTS(DUMMY%J))%DATA%X(6*(K1-1)+2)
                 C(3,K1)=TS(MI(DUMMY%I)%DATA%CTS(DUMMY%J))%DATA%X(6*(K1-1)+3)
                 C(1,K1+NRES)=TS(MI(DUMMY%I)%DATA%CTS(DUMMY%J))%DATA%X(6*(K1-1)+4)
                 C(2,K1+NRES)=TS(MI(DUMMY%I)%DATA%CTS(DUMMY%J))%DATA%X(6*(K1-1)+5)
                 C(3,K1+NRES)=TS(MI(DUMMY%I)%DATA%CTS(DUMMY%J))%DATA%X(6*(K1-1)+6)
              ENDDO
              CALL UPDATEDC
              CALL INT_FROM_CART(.TRUE.,.FALSE.)
              CALL CHAINBUILD
           ENDIF
!          CALL POTENTIAL(TS(MI(DUMMY%I)%DATA%CTS(DUMMY%J))%DATA%X,DUMMY1,GRAD,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
           IF (CHRMMT) THEN
              HORDER=1
              FPGRP='C1'
              IF (MACHINE) THEN
                 WRITE(88) HORDER,FPGRP
              ELSE
                 WRITE(88,'(I6,1X,A4)') HORDER,FPGRP
              ENDIF
              IF (.NOT.NOFRQS) THEN
                 IF (ENDNUMHESS) THEN
                    CALL MAKENUMHESS(TS(MI(DUMMY%I)%DATA%CTS(DUMMY%J))%DATA%X,NATOMS)
                 ELSE
                    CALL POTENTIAL(TS(MI(DUMMY%I)%DATA%CTS(DUMMY%J))%DATA%X,DUMMY1,GRAD,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
                 ENDIF
                 CALL MASSWT2(NATOMS,ATMASS,DUMQ,GRAD,.TRUE.)
                 CALL DSYEV('N','U',3*NATOMS,HESS,SIZE(HESS,1),DIAG,TEMPA,9*NATOMS,INFO)
                 IF (DIAG(1).LT.DIAG(3*NATOMS)) CALL EIGENSORT_VAL_ASC(DIAG,HESS,3*NATOMS,3*NATOMS)
                 IF (MACHINE) THEN
                    WRITE(88) (DIAG(J2)*4.184D26,J2=1,3*NATOMS)
                 ELSE
                    WRITE(88,'(3G20.10)') (DIAG(J2)*4.184D26,J2=1,3*NATOMS)
                 ENDIF
              ENDIF
           ELSEIF (UNRST) THEN
              HORDER=1
              FPGRP='C1'
              WRITE(88,'(I6,1X,A4)') HORDER,FPGRP
              IF (.NOT.NOFRQS) THEN
                 IF (ENDNUMHESS) THEN
                    CALL MAKENUMINTHESS(NINTS,NATOMS)
                    CALL GETSTUFF(KD,NNZ,NINTB)
                    CALL INTSECDET(TS(MI(DUMMY%I)%DATA%CTS(DUMMY%J))%DATA%X,3*NATOMS,KD,NNZ,NINTB,DIAG)
                 ELSE
                    CALL POTENTIAL(TS(MI(DUMMY%I)%DATA%CTS(DUMMY%J))%DATA%X,DUMMY1,GRAD,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
                 ENDIF
                 DO J2=1,NINTS-1
                    IF (DIAG(J2).LT.0.0D0) PRINT *,'HIGHER ORDER SADDLE FOUND IN PATHWAY - TS ',I,'EIGENVALUE ',DIAG(J2)
                 END DO
                 WRITE(88,'(3G20.10)') (DIAG(J2),J2=1,3*NATOMS)
              ENDIF
           ELSEIF (AMHT) THEN
              PRINT '(A)','IN TS AMH BLOCK WRITING 1 C1'
              WRITE(88,'(I6,1X,A4)') 1,'C1'
              IF (.NOT.NOFRQS) THEN
              IF (ENDNUMHESS) THEN
                CALL MAKENUMINTHESS(NINTS,NATOMS)
                CALL GETSTUFF(KD,NNZ,NINTB)
                CALL INTSECDET(TS(MI(DUMMY%I)%DATA%CTS(DUMMY%J))%DATA%X,3*NATOMS,KD,NNZ,NINTB,DIAG)
             ELSE
                CALL POTENTIAL(TS(MI(DUMMY%I)%DATA%CTS(DUMMY%J))%DATA%X,DUMMY1,GRAD,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
              ENDIF
              DO J2=1,NINTS-1
                 IF (DIAG(J2).LT.0.0D0) PRINT *,'HIGHER ORDER SADDLE FOUND IN PATHWAY - TS ',I,'EIGENVALUE ',DIAG(J2)
              END DO
              WRITE(88,'(3G20.10)') (DIAG(J2),J2=1,3*NATOMS)
           ENDIF

        ELSE
           CALL SYMMETRY(HORDER,.FALSE.,TS(MI(DUMMY%I)%DATA%CTS(DUMMY%J))%DATA%X,INERTIA)
              WRITE(88,'(I6,1X,A4)') HORDER,FPGRP
              IF (.NOT.NOFRQS) THEN
                 IF (ENDNUMHESS) THEN
                    CALL MAKENUMHESS(TS(MI(DUMMY%I)%DATA%CTS(DUMMY%J))%DATA%X,NATOMS)
                 ELSE
                    CALL POTENTIAL(TS(MI(DUMMY%I)%DATA%CTS(DUMMY%J))%DATA%X,DUMMY1,GRAD,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
                 ENDIF
                 CALL MASSWT(NATOMS,ATMASS,DUMQ,GRAD,.TRUE.)
                 CALL DSYEV('N','U',3*NATOMS,HESS,SIZE(HESS,1),DIAG,TEMPA,9*NATOMS,INFO)
                 IF (DIAG(1).LT.DIAG(3*NATOMS)) CALL EIGENSORT_VAL_ASC(DIAG,HESS,3*NATOMS,3*NATOMS)
                 IF (SDT) THEN
                    WRITE(88,'(3G20.10)') (DIAG(J2)*4.184D26,J2=1,3*NATOMS)
                 ELSE
                    WRITE(88,'(3G20.10)') (DIAG(J2),J2=1,3*NATOMS)
                 ENDIF
              ENDIF
           ENDIF
        ELSE
           CALL SYMMETRY(HORDER,.FALSE.,TS(MI(DUMMY%I)%DATA%CTS(DUMMY%J))%DATA%X,INERTIA)
           WRITE(88,'(I6,1X,A4)') HORDER,FPGRP
!          WRITE(88,'(3G20.10)') (FRQSTS(J2,WHICHTS(J1)),J2=1,3*NATOMS)
        ENDIF
        CALL MINPERMDIST(MI(DUMMY%I)%DATA%X,TS(MI(DUMMY%I)%DATA%CTS(DUMMY%J))%DATA%X,NATOMS,DEBUG, &
    &                    PARAM1,PARAM2,PARAM3,BULKT,TWOD,DIST,DIST2,RIGIDBODY,RMAT)
        DISTPF=DIST
        IF (MACHINE) THEN
           WRITE(88) (TS(MI(DUMMY%I)%DATA%CTS(DUMMY%J))%DATA%X(J2),J2=1,NOPT)
        ELSE
           WRITE(88,'(3F25.15)') (TS(MI(DUMMY%I)%DATA%CTS(DUMMY%J))%DATA%X(J2),J2=1,NOPT)
        ENDIF
        PREVIOUSTS=TS(MI(DUMMY%I)%DATA%CTS(DUMMY%J))%DATA%X
        DUMMY=>DUMMY%NEXT
        I=I+1
     ENDDO
     CALL FLUSH(88,ISTAT)
     CLOSE(88)

     END SUBROUTINE MAKEPATHINFO

!
!  DUMP THE LATEST MIN-SAD-MIN TRIPLE TO PATH.INFO IN THE USUAL FORMAT
!  
     SUBROUTINE MAKEALLPATHINFO(QTS,QPLUS,QMINUS,ETS,EPLUS,EMINUS,FRQSTS,FRQSPLUS,FRQSMINUS)
     USE SYMINF 
     USE MODHESS
     USE MODCHARMM
     USE MODUNRES
     USE KEY, ONLY: FILTH, FILTHSTR, UNRST, TWOD, BULKT, MACHINE, NOFRQS, AMBERT, NABT, AMHT, SEQ, TARFL, NRES_AMH_TEMP, SDT
     USE COMMONS, ONLY: ATMASS, NINTS, ZSYM
     USE PORFUNCS
     IMPLICIT NONE

     CHARACTER(LEN=20) :: PINFOSTRING
     DOUBLE PRECISION :: DIHE,ALLANG,DISTPF,DUMMY1,GRAD(3*NATOMS),RMS,DIAG(3*NATOMS),TEMPA(9*NATOMS),DUMQ(3*NATOMS)
     INTEGER :: HORDER,I,INFO,J2,K1,RECLEN,ISTAT,J1
     LOGICAL :: BTEST,KD,NNZ,NINTB,TSFRQDONE,MINFRQDONE
     DOUBLE PRECISION :: QTS(3*NATOMS),QPLUS(3*NATOMS),QMINUS(3*NATOMS),FRQSTS(3*NATOMS),FRQSPLUS(3*NATOMS),FRQSMINUS(3*NATOMS), &
    &                    ETS,EPLUS,EMINUS,INERTIA(3,3)

!    LOCAL AMH VARIABLES
     INTEGER :: I_RES, GLY_COUNT

      LOGICAL KNOWE, KNOWG, KNOWH
      COMMON /KNOWN/ KNOWE, KNOWG, KNOWH
!   
!  FILE PATH.INFO IS NOW OPENED IN KEYWORD.F
!  
!    IF (FILTH.EQ.0) THEN
!       WRITE(PINFOSTRING,'(A9)') 'PATH.INFO'
!    ELSE
!       WRITE(PINFOSTRING,'(A)') 'PATH.INFO.'//TRIM(ADJUSTL(FILTHSTR))
!    ENDIF
!    IF (MACHINE) THEN
!         OPEN(UNIT=88,FILE=PINFOSTRING,STATUS='UNKNOWN',FORM='UNFORMATTED',POSITION='APPEND')
!    ELSE
!         OPEN(UNIT=88,FILE=PINFOSTRING,STATUS='UNKNOWN',POSITION='APPEND')
!    ENDIF

     TSFRQDONE=.FALSE.  ! ASSUME THAT WE NEVER KNOW THE FREQUENCIES
     MINFRQDONE=.FALSE. ! ASSUME THAT WE NEVER KNOW THE FREQUENCIES
!
!  FIRST DUMP THE + MINIMUM.
!  
     IF (UNRST.AND.CALCDIHE) THEN
        CALL UNRESCALCDIHEREF(DIHE,ALLANG,QPLUS)
        WRITE(88,'(3F20.10)') EPLUS, DIHE
     ELSE
        WRITE(88,'(F20.10)') EPLUS
     ENDIF
     IF ((.NOT.MINFRQDONE).OR.(TWOD).OR.CHRMMT.OR.UNRST.OR.AMBERT.OR.NABT) THEN
        IF (UNRST) THEN ! JMC UPDATE COORDS
           DO K1=1,NRES 
               C(1,K1)=QPLUS(6*(K1-1)+1)
               C(2,K1)=QPLUS(6*(K1-1)+2)
               C(3,K1)=QPLUS(6*(K1-1)+3)
               C(1,K1+NRES)=QPLUS(6*(K1-1)+4)
               C(2,K1+NRES)=QPLUS(6*(K1-1)+5)
               C(3,K1+NRES)=QPLUS(6*(K1-1)+6)
            ENDDO
            CALL UPDATEDC
            CALL INT_FROM_CART(.TRUE.,.FALSE.)
            CALL CHAINBUILD
        ENDIF
        IF (CHRMMT) NCHENCALLS = 999
!       CALL POTENTIAL(QPLUS,DUMMY1,GRAD,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
        IF (CHRMMT) THEN
           HORDER=1 
           FPGRP='C1'
           IF (MACHINE) THEN
              WRITE(88) HORDER,FPGRP
           ELSE
              WRITE(88,'(I6,1X,A4)') HORDER,FPGRP
           ENDIF
           IF (.NOT.NOFRQS) THEN
              IF (ENDNUMHESS) THEN
                 CALL MAKENUMHESS(QPLUS,NATOMS)
              ELSE
                 CALL POTENTIAL(QPLUS,DUMMY1,GRAD,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
              ENDIF
              CALL MASSWT2(NATOMS,ATMASS,DUMQ,GRAD,.TRUE.)
              CALL DSYEV('N','U',3*NATOMS,HESS,SIZE(HESS,1),DIAG,TEMPA,9*NATOMS,INFO)
              IF (DIAG(1).LT.DIAG(3*NATOMS)) CALL EIGENSORT_VAL_ASC(DIAG,HESS,3*NATOMS,3*NATOMS)
              IF (MACHINE) THEN
                 WRITE(88) (DIAG(J2)*4.184D26,J2=1,3*NATOMS)
              ELSE
                 WRITE(88,'(3G20.10)') (DIAG(J2)*4.184D26,J2=1,3*NATOMS)
              ENDIF
           ENDIF
        ELSEIF (AMBERT.OR.NABT) THEN
           HORDER=1 
           FPGRP='C1'
           IF (MACHINE) THEN
              WRITE(88) HORDER,FPGRP
           ELSE
              WRITE(88,'(I6,1X,A4)') HORDER,FPGRP
           ENDIF
           IF (.NOT.NOFRQS) THEN
              IF (ENDNUMHESS.OR.AMBERT) THEN
                 CALL MAKENUMHESS(QPLUS,NATOMS)
              ELSE
                 CALL POTENTIAL(QPLUS,DUMMY1,GRAD,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
              ENDIF
              CALL MASSWT2(NATOMS,ATMASS,DUMQ,GRAD,.TRUE.)
              CALL DSYEV('N','U',3*NATOMS,HESS,SIZE(HESS,1),DIAG,TEMPA,9*NATOMS,INFO)
              IF (DIAG(1).LT.DIAG(3*NATOMS)) CALL EIGENSORT_VAL_ASC(DIAG,HESS,3*NATOMS,3*NATOMS)
              IF (MACHINE) THEN
                 WRITE(88) (DIAG(J2)*4.184D26,J2=1,3*NATOMS)
              ELSE
                 WRITE(88,'(3G20.10)') (DIAG(J2)*4.184D26,J2=1,3*NATOMS)
              ENDIF
           ENDIF
        ELSEIF (UNRST) THEN
           HORDER=1
           FPGRP='C1'
           WRITE(88,'(I6,1X,A4)') HORDER,FPGRP
           IF (.NOT.NOFRQS) THEN
              IF (ENDNUMHESS) THEN
                 CALL MAKENUMINTHESS(NINTS,NATOMS)
                 CALL GETSTUFF(KD,NNZ,NINTB)
                 CALL INTSECDET(QPLUS,3*NATOMS,KD,NNZ,NINTB,DIAG)
              ELSE
                 CALL POTENTIAL(QPLUS,DUMMY1,GRAD,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
              ENDIF
              WRITE(88,'(3G20.10)') (DIAG(J2),J2=1,3*NATOMS)
           ENDIF
        ELSEIF (AMHT) THEN
           WRITE(88,'(I6,1X,A4)') 1,' C1'
           IF (.NOT.NOFRQS) THEN
              IF (ENDNUMHESS) THEN
                 CALL MAKENUMHESS(QPLUS,NATOMS)
              ELSE
                 CALL POTENTIAL(QPLUS,DUMMY1,GRAD,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
              ENDIF
            CALL MASSWT(NATOMS,ATMASS,DUMQ,GRAD,.TRUE.)
            CALL DSYEV('N','U',3*NATOMS,HESS,SIZE(HESS,1),DIAG,TEMPA,9*NATOMS,INFO)
            IF (DIAG(1).LT.DIAG(3*NATOMS)) CALL EIGENSORT_VAL_ASC(DIAG,HESS,3*NATOMS,3*NATOMS)
             WRITE(88,'(3G20.10)') (DIAG(J2),J2=1,3*NATOMS)
         ENDIF
        ELSE
           CALL SYMMETRY(HORDER,.FALSE.,QPLUS,INERTIA)
           WRITE(88,'(I6,1X,A4)') HORDER,FPGRP
           IF (.NOT.NOFRQS) THEN
              IF (ENDNUMHESS) THEN
                 CALL MAKENUMHESS(QPLUS,NATOMS)
              ELSE
                 CALL POTENTIAL(QPLUS,DUMMY1,GRAD,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
              ENDIF
              CALL MASSWT(NATOMS,ATMASS,DUMQ,GRAD,.TRUE.)
              CALL DSYEV('N','U',3*NATOMS,HESS,SIZE(HESS,1),DIAG,TEMPA,9*NATOMS,INFO)
              IF (DIAG(1).LT.DIAG(3*NATOMS)) CALL EIGENSORT_VAL_ASC(DIAG,HESS,3*NATOMS,3*NATOMS)
              IF (SDT) THEN
                 WRITE(88,'(3G20.10)') (DIAG(J2)*4.184D26,J2=1,3*NATOMS)
              ELSE
                 WRITE(88,'(3G20.10)') (DIAG(J2),J2=1,3*NATOMS)
              ENDIF
           ENDIF
        ENDIF
     ELSE
        CALL SYMMETRY(HORDER,.FALSE.,QPLUS,INERTIA)
        WRITE(88,'(I6,1X,A4)') HORDER,FPGRP
! JMC   WRITE(88,'(3G20.10)') (FSAVEMIN(J2,J1),J2=1,3*NATOMS)
     ENDIF
     IF (MACHINE) THEN
          WRITE(88) (QPLUS,J2=1,3*NATOMS)
     ELSEIF (AMHT) THEN

!  THIS IS FOR PLACE HOLDING C-BETAS FOR GLYCINE IN AMH
        GLY_COUNT = 0

        DO J2=1, NRES_AMH_TEMP
           IF (SEQ(J2).EQ.8) THEN
!             WRITE(2,*)SEQ(J2) , J2
               WRITE(88,*)QPLUS(9*(J2-1)+1-GLY_COUNT*3), &
               QPLUS(9*(J2-1)+2-GLY_COUNT*3),QPLUS(9*(J2-1)+3-GLY_COUNT*3)
              WRITE(88,*)QPLUS(9*(J2-1)+1-GLY_COUNT*3), &
               QPLUS(9*(J2-1)+2-GLY_COUNT*3),QPLUS(9*(J2-1)+3-GLY_COUNT*3)
              WRITE(88,*)QPLUS(9*(J2-1)+4-GLY_COUNT*3), &
                QPLUS(9*(J2-1)+5-GLY_COUNT*3),QPLUS(9*(J2-1)+6-GLY_COUNT*3)
              GLY_COUNT = GLY_COUNT + 1
           ELSE
!            WRITE(2,*)SEQ(J2) , J2
             WRITE(88,*)QPLUS(9*(J2-1)+1-GLY_COUNT*3), &
               QPLUS(9*(J2-1)+2-GLY_COUNT*3),QPLUS(9*(J2-1)+3-GLY_COUNT*3)
             WRITE(88,*)QPLUS(9*(J2-1)+4-GLY_COUNT*3), &
               QPLUS(9*(J2-1)+5-GLY_COUNT*3),QPLUS(9*(J2-1)+6-GLY_COUNT*3)
            WRITE(88,*)QPLUS(9*(J2-1)+7-GLY_COUNT*3), &
               QPLUS(9*(J2-1)+8-GLY_COUNT*3),QPLUS(9*(J2-1)+9-GLY_COUNT*3)
           ENDIF
        ENDDO
      ELSE
          WRITE(88,'(3F25.15)') (QPLUS(J2),J2=1,3*NATOMS)
     ENDIF
!
! NOW THE TRANSITION STATE
!
     IF (MACHINE) THEN
          WRITE(88) ETS
     ELSE
          WRITE(88,'(F20.10)') ETS
     ENDIF
     IF ((.NOT.TSFRQDONE).OR.(TWOD).OR.CHRMMT.OR.UNRST.OR.AMBERT.OR.NABT) THEN
        IF (UNRST) THEN
          DO K1=1,NRES ! JMC UPDATE COORDS
              C(1,K1)=QTS(6*(K1-1)+1)
              C(2,K1)=QTS(6*(K1-1)+2)
              C(3,K1)=QTS(6*(K1-1)+3)
              C(1,K1+NRES)=QTS(6*(K1-1)+4)
              C(2,K1+NRES)=QTS(6*(K1-1)+5)
              C(3,K1+NRES)=QTS(6*(K1-1)+6)
           ENDDO
           CALL UPDATEDC
           CALL INT_FROM_CART(.TRUE.,.FALSE.)
           CALL CHAINBUILD
        ENDIF
        IF (CHRMMT) NCHENCALLS = 999
!       CALL POTENTIAL(QTS,DUMMY1,GRAD,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
        IF (CHRMMT) THEN
           HORDER=1
           FPGRP='C1'
           IF (MACHINE) THEN
              WRITE(88) HORDER,FPGRP
           ELSE
              WRITE(88,'(I6,1X,A4)') HORDER,FPGRP
           ENDIF
           IF (.NOT.NOFRQS) THEN
              IF (ENDNUMHESS) THEN
                 CALL MAKENUMHESS(QTS,NATOMS)
              ELSE
                 CALL POTENTIAL(QTS,DUMMY1,GRAD,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
              ENDIF
              CALL MASSWT2(NATOMS,ATMASS,DUMQ,GRAD,.TRUE.)
              CALL DSYEV('N','U',3*NATOMS,HESS,SIZE(HESS,1),DIAG,TEMPA,9*NATOMS,INFO)
              IF (DIAG(1).LT.DIAG(3*NATOMS)) CALL EIGENSORT_VAL_ASC(DIAG,HESS,3*NATOMS,3*NATOMS)
              IF (MACHINE) THEN
                 WRITE(88) (DIAG(J2)*4.184D26,J2=1,3*NATOMS)
              ELSE
                 WRITE(88,'(3G20.10)') (DIAG(J2)*4.184D26,J2=1,3*NATOMS)
              ENDIF
           ENDIF
        ELSE IF (AMBERT.OR.NABT) THEN
           HORDER=1
           FPGRP='C1'
           IF (MACHINE) THEN
              WRITE(88) HORDER,FPGRP
           ELSE
              WRITE(88,'(I6,1X,A4)') HORDER,FPGRP
           ENDIF
           IF (.NOT.NOFRQS) THEN
              IF (ENDNUMHESS.OR.AMBERT) THEN
                 CALL MAKENUMHESS(QTS,NATOMS)
              ELSE
                 CALL POTENTIAL(QTS,DUMMY1,GRAD,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
              ENDIF
              CALL MASSWT2(NATOMS,ATMASS,DUMQ,GRAD,.TRUE.)
!              OPEN(UNIT=955, FILE='HESSDUMP', STATUS='UNKNOWN')
!              WRITE(955,*) HESS(:,:)
              CALL DSYEV('N','U',3*NATOMS,HESS,SIZE(HESS,1),DIAG,TEMPA,9*NATOMS,INFO)
              IF (DIAG(1).LT.DIAG(3*NATOMS)) CALL EIGENSORT_VAL_ASC(DIAG,HESS,3*NATOMS,3*NATOMS)
               
              IF (MACHINE) THEN
                 WRITE(88) (DIAG(J2)*4.184D26,J2=1,3*NATOMS)
              ELSE
                 WRITE(88,'(3G20.10)') (DIAG(J2)*4.184D26,J2=1,3*NATOMS)
              ENDIF
           ENDIF
        ELSEIF (UNRST) THEN
           HORDER=1
           FPGRP='C1'
           WRITE(88,'(I6,1X,A4)') HORDER,FPGRP
           IF (.NOT.NOFRQS) THEN
              IF (ENDNUMHESS) THEN
                 CALL MAKENUMINTHESS(NINTS,NATOMS)
                 CALL GETSTUFF(KD,NNZ,NINTB)
                 CALL INTSECDET(QTS,3*NATOMS,KD,NNZ,NINTB,DIAG)
              ELSE
                 CALL POTENTIAL(QTS,DUMMY1,GRAD,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
              ENDIF
              DO J2=1,NINTS-1
                 IF (DIAG(J2).LT.0.0D0) PRINT *,'HIGHER ORDER SADDLE FOUND IN PATHWAY - TS ',I,'EIGENVALUE ',DIAG(J2)
              END DO
              WRITE(88,'(3G20.10)') (DIAG(J2),J2=1,3*NATOMS)
           ENDIF
        ELSEIF (AMHT) THEN
           WRITE(88,'(I6,1X,A4)') 1,' C1'
           IF (.NOT.NOFRQS) THEN
              IF (ENDNUMHESS) THEN
                 CALL MAKENUMHESS(QTS,NATOMS)
              ELSE
                 CALL POTENTIAL(QTS,DUMMY1,GRAD,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
              ENDIF
              CALL MASSWT(NATOMS,ATMASS,DUMQ,GRAD,.TRUE.)
              CALL DSYEV('N','U',3*NATOMS,HESS,SIZE(HESS,1),DIAG,TEMPA,9*NATOMS,INFO)
              IF (DIAG(1).LT.DIAG(3*NATOMS)) CALL EIGENSORT_VAL_ASC(DIAG,HESS,3*NATOMS,3*NATOMS)
              WRITE(88,'(3G20.10)') (DIAG(J2),J2=1,3*NATOMS)
           ENDIF
        ELSE
           CALL SYMMETRY(HORDER,.FALSE.,QTS,INERTIA)
           WRITE(88,'(I6,1X,A4)') HORDER,FPGRP
           IF (.NOT.NOFRQS) THEN
              IF (ENDNUMHESS) THEN
                 CALL MAKENUMHESS(QTS,NATOMS)
              ELSE
                 CALL POTENTIAL(QTS,DUMMY1,GRAD,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
              ENDIF
              CALL MASSWT(NATOMS,ATMASS,DUMQ,GRAD,.TRUE.)
              CALL DSYEV('N','U',3*NATOMS,HESS,SIZE(HESS,1),DIAG,TEMPA,9*NATOMS,INFO)
              IF (DIAG(1).LT.DIAG(3*NATOMS)) CALL EIGENSORT_VAL_ASC(DIAG,HESS,3*NATOMS,3*NATOMS)
              IF (SDT) THEN
                 WRITE(88,'(3G20.10)') (DIAG(J2)*4.184D26,J2=1,3*NATOMS)
              ELSE
                 WRITE(88,'(3G20.10)') (DIAG(J2),J2=1,3*NATOMS)
              ENDIF
           ENDIF
        ENDIF
     ELSE
        CALL SYMMETRY(HORDER,.FALSE.,QTS,INERTIA)
        WRITE(88,'(I6,1X,A4)') HORDER,FPGRP
!             WRITE(88,'(3G20.10)') (FRQSTS(J2,WHICHTS(J1)),J2=1,3*NATOMS)
     ENDIF
     IF (MACHINE) THEN
          WRITE(88) (QTS(J2),J2=1,NOPT)
     ELSEIF (AMHT) THEN
!       READ SEQUENCE

!  THIS IS FOR PLACE HOLDING C-BETAS FOR GLYCINE IN AMH
        GLY_COUNT = 0

        DO J2=1, NRES_AMH_TEMP
           IF (SEQ(J2).EQ.8) THEN
!             WRITE(2,*)SEQ(J2) , J2
               WRITE(88,*)QTS(9*(J2-1)+1-GLY_COUNT*3), &
                QTS(9*(J2-1)+2-GLY_COUNT*3),QTS(9*(J2-1)+3-GLY_COUNT*3)
              WRITE(88,*)QTS(9*(J2-1)+1-GLY_COUNT*3), &
                QTS(9*(J2-1)+2-GLY_COUNT*3),QTS(9*(J2-1)+3-GLY_COUNT*3)
              WRITE(88,*)QTS(9*(J2-1)+4-GLY_COUNT*3), &
                QTS(9*(J2-1)+5-GLY_COUNT*3),QTS(9*(J2-1)+6-GLY_COUNT*3)
              GLY_COUNT = GLY_COUNT + 1
           ELSE
!            WRITE(2,*)SEQ(J2) , J2
             WRITE(88,*)QTS(9*(J2-1)+1-GLY_COUNT*3), &
               QTS(9*(J2-1)+2-GLY_COUNT*3),QTS(9*(J2-1)+3-GLY_COUNT*3)
             WRITE(88,*)QTS(9*(J2-1)+4-GLY_COUNT*3), &
               QTS(9*(J2-1)+5-GLY_COUNT*3),QTS(9*(J2-1)+6-GLY_COUNT*3)
             WRITE(88,*)QTS(9*(J2-1)+7-GLY_COUNT*3), &
               QTS(9*(J2-1)+8-GLY_COUNT*3),QTS(9*(J2-1)+9-GLY_COUNT*3)
           ENDIF
         ENDDO
     ELSE
          WRITE(88,'(3F25.15)') (QTS(J2),J2=1,NOPT)
     ENDIF
!
!  FINALLY DUMP THE - MINIMUM.
!
     IF (UNRST.AND.CALCDIHE) THEN
        CALL UNRESCALCDIHEREF(DIHE,ALLANG,QMINUS)
        WRITE(88,'(3F20.10)') EMINUS, DIHE
     ELSE
        WRITE(88,'(F20.10)') EMINUS
     ENDIF
     IF ((.NOT.MINFRQDONE).OR.(TWOD).OR.CHRMMT.OR.UNRST.OR.AMBERT.OR.NABT) THEN
        IF (UNRST) THEN ! JMC UPDATE COORDS
           DO K1=1,NRES
               C(1,K1)=QMINUS(6*(K1-1)+1)
               C(2,K1)=QMINUS(6*(K1-1)+2)
               C(3,K1)=QMINUS(6*(K1-1)+3)
               C(1,K1+NRES)=QMINUS(6*(K1-1)+4)
               C(2,K1+NRES)=QMINUS(6*(K1-1)+5)
               C(3,K1+NRES)=QMINUS(6*(K1-1)+6)
            ENDDO
            CALL UPDATEDC
            CALL INT_FROM_CART(.TRUE.,.FALSE.)
            CALL CHAINBUILD
        ENDIF
        IF (CHRMMT) NCHENCALLS = 999
!       CALL POTENTIAL(QMINUS,DUMMY1,GRAD,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
        IF (CHRMMT) THEN
           HORDER=1 
           FPGRP='C1'
           IF (MACHINE) THEN
              WRITE(88) HORDER,FPGRP
           ELSE
              WRITE(88,'(I6,1X,A4)') HORDER,FPGRP
           ENDIF
           IF (.NOT.NOFRQS) THEN
              IF (ENDNUMHESS) THEN
                 CALL MAKENUMHESS(QMINUS,NATOMS)
              ELSE
                 CALL POTENTIAL(QMINUS,DUMMY1,GRAD,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
              ENDIF
              CALL MASSWT2(NATOMS,ATMASS,DUMQ,GRAD,.TRUE.)
              CALL DSYEV('N','U',3*NATOMS,HESS,SIZE(HESS,1),DIAG,TEMPA,9*NATOMS,INFO)
              IF (DIAG(1).LT.DIAG(3*NATOMS)) CALL EIGENSORT_VAL_ASC(DIAG,HESS,3*NATOMS,3*NATOMS)
              IF (MACHINE) THEN
                   WRITE(88) (DIAG(J2)*4.184D26,J2=1,3*NATOMS)
              ELSE
                   WRITE(88,'(3G20.10)') (DIAG(J2)*4.184D26,J2=1,3*NATOMS)
              ENDIF
           ENDIF
        ELSEIF (AMBERT.OR.NABT) THEN
           HORDER=1
           FPGRP='C1'
           IF (MACHINE) THEN
              WRITE(88) HORDER,FPGRP
           ELSE
              WRITE(88,'(I6,1X,A4)') HORDER,FPGRP
           ENDIF
           IF (.NOT.NOFRQS) THEN
              IF (ENDNUMHESS.OR.AMBERT) THEN
                 CALL MAKENUMHESS(QMINUS,NATOMS)
              ELSE
                 CALL POTENTIAL(QMINUS,DUMMY1,GRAD,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
              ENDIF
              CALL MASSWT2(NATOMS,ATMASS,DUMQ,GRAD,.TRUE.)
              CALL DSYEV('N','U',3*NATOMS,HESS,SIZE(HESS,1),DIAG,TEMPA,9*NATOMS,INFO)
              IF (DIAG(1).LT.DIAG(3*NATOMS)) CALL EIGENSORT_VAL_ASC(DIAG,HESS,3*NATOMS,3*NATOMS)
              IF (MACHINE) THEN
                 WRITE(88) (DIAG(J2)*4.184D26,J2=1,3*NATOMS)
              ELSE
                 WRITE(88,'(3G20.10)') (DIAG(J2)*4.184D26,J2=1,3*NATOMS)
              ENDIF
           ENDIF
        ELSEIF (UNRST) THEN
           HORDER=1
           FPGRP='C1'
           WRITE(88,'(I6,1X,A4)') HORDER,FPGRP
           IF (.NOT.NOFRQS) THEN
              IF (ENDNUMHESS) THEN
                 CALL MAKENUMINTHESS(NINTS,NATOMS)
                 CALL GETSTUFF(KD,NNZ,NINTB)
                 CALL INTSECDET(QMINUS,3*NATOMS,KD,NNZ,NINTB,DIAG)
              ELSE
                 CALL POTENTIAL(QMINUS,DUMMY1,GRAD,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
              ENDIF
              WRITE(88,'(3G20.10)') (DIAG(J2),J2=1,3*NATOMS)
           ENDIF
        ELSEIF (AMHT) THEN
           WRITE(88,'(I6,1X,A4)') 1,' C1'
           IF (.NOT.NOFRQS) THEN
              IF (ENDNUMHESS) THEN
                 CALL MAKENUMHESS(QMINUS,NATOMS)
              ELSE
                 CALL POTENTIAL(QMINUS,DUMMY1,GRAD,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
             ENDIF
              CALL MASSWT(NATOMS,ATMASS,DUMQ,GRAD,.TRUE.)
              CALL DSYEV('N','U',3*NATOMS,HESS,SIZE(HESS,1),DIAG,TEMPA,9*NATOMS,INFO)
              IF (DIAG(1).LT.DIAG(3*NATOMS)) CALL EIGENSORT_VAL_ASC(DIAG,HESS,3*NATOMS,3*NATOMS)
              WRITE(88,'(3G20.10)') (DIAG(J2),J2=1,3*NATOMS)
           ENDIF
        ELSE
           CALL SYMMETRY(HORDER,.FALSE.,QMINUS,INERTIA)
           WRITE(88,'(I6,1X,A4)') HORDER,FPGRP
           IF (.NOT.NOFRQS) THEN
              IF (ENDNUMHESS) THEN
                 CALL MAKENUMHESS(QMINUS,NATOMS)
              ELSE
                 CALL POTENTIAL(QMINUS,DUMMY1,GRAD,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
              ENDIF
              CALL MASSWT(NATOMS,ATMASS,DUMQ,GRAD,.TRUE.)
              CALL DSYEV('N','U',3*NATOMS,HESS,SIZE(HESS,1),DIAG,TEMPA,9*NATOMS,INFO)
              IF (DIAG(1).LT.DIAG(3*NATOMS)) CALL EIGENSORT_VAL_ASC(DIAG,HESS,3*NATOMS,3*NATOMS)
              IF (SDT) THEN
                 WRITE(88,'(3G20.10)') (DIAG(J2)*4.184D26,J2=1,3*NATOMS)
              ELSE
                 WRITE(88,'(3G20.10)') (DIAG(J2),J2=1,3*NATOMS)
              ENDIF
           ENDIF
        ENDIF
     ELSE
        CALL SYMMETRY(HORDER,.FALSE.,QMINUS,INERTIA)
        WRITE(88,'(I6,1X,A4)') HORDER,FPGRP
! JMC   WRITE(88,'(3G20.10)') (FSAVEMIN(J2,J1),J2=1,3*NATOMS)
     ENDIF
     IF (MACHINE) THEN
        WRITE(88) (QMINUS,J2=1,3*NATOMS)
     ELSEIF (AMHT) THEN
!       READ SEQUENCE

!  THIS IS FOR PLACE HOLDING C-BETAS FOR GLYCINE IN AMH
        GLY_COUNT = 0

        DO J2=1,NRES_AMH_TEMP
           IF (SEQ(J2).EQ.8) THEN
!             WRITE(2,*)SEQ(J2) , J2
               WRITE(88,*)QMINUS(9*(J2-1)+1-GLY_COUNT*3), &
                QMINUS(9*(J2-1)+2-GLY_COUNT*3),QMINUS(9*(J2-1)+3-GLY_COUNT*3)
              WRITE(88,*)QMINUS(9*(J2-1)+1-GLY_COUNT*3), &
                QMINUS(9*(J2-1)+2-GLY_COUNT*3),QMINUS(9*(J2-1)+3-GLY_COUNT*3)
              WRITE(88,*)QMINUS(9*(J2-1)+4-GLY_COUNT*3), &
               QMINUS(9*(J2-1)+5-GLY_COUNT*3),QMINUS(9*(J2-1)+6-GLY_COUNT*3)
              GLY_COUNT = GLY_COUNT + 1
           ELSE
!            WRITE(2,*)SEQ(J2) , J2
             WRITE(88,*)QMINUS(9*(J2-1)+1-GLY_COUNT*3), &
               QMINUS(9*(J2-1)+2-GLY_COUNT*3),QMINUS(9*(J2-1)+3-GLY_COUNT*3)
             WRITE(88,*)QMINUS(9*(J2-1)+4-GLY_COUNT*3), &
               QMINUS(9*(J2-1)+5-GLY_COUNT*3),QMINUS(9*(J2-1)+6-GLY_COUNT*3)
            WRITE(88,*)QMINUS(9*(J2-1)+7-GLY_COUNT*3), &
               QMINUS(9*(J2-1)+8-GLY_COUNT*3),QMINUS(9*(J2-1)+9-GLY_COUNT*3)
           ENDIF
       ENDDO
     ELSE
        WRITE(88,'(3F25.15)') (QMINUS(J2),J2=1,3*NATOMS)
     ENDIF

     KNOWH = .FALSE. ! NEEDED OTHERWISE THE NEXT TS SEARCH WILL USE THE WRONG HESSIAN, IF ONE IS REQUIRED.
     CALL FLUSH(88,ISTAT)

     END SUBROUTINE MAKEALLPATHINFO

     SUBROUTINE CHECKPAIR(I,J,PERMTEST) ! CHECKS THAT MINIMA I AND J ARE DIFFERENT AND PUTS THEM INTO THE ORIENTATION WITH MINIMAL D
          USE PORFUNCS
          USE KEYUTILS
          USE KEY,ONLY : DEBUG, RIGIDBODY, PERMDIST,TWOD,BULKT,NEBRESEEDT,LOCALPERMDIST
          USE COMMONS,ONLY : PARAM1,PARAM2,PARAM3
          IMPLICIT NONE
          DOUBLE PRECISION RMAT(3,3), DIST2, DIST
          INTEGER,INTENT(IN) :: I,J
          INTEGER J1, J2
          LOGICAL PERMTEST
    
          PERMTEST=.FALSE.
          IF (ABS(MI(I)%DATA%E-MI(J)%DATA%E) < EDIFFTOL) THEN
             PERMTEST=.TRUE.
             WRITE(*,'(/1X,A,2I5,A)') "ENERGIES OF THE MINIMA IN THE PAIR ",I,J," ARE THE SAME - CHECKING DISTANCE ..."

!            IF (BULKT) THEN
!               CALL NEWMINDIST(MI(I)%DATA%X,MI(J)%DATA%X,NATOMS,D,BULKT,TWOD,'AX   ',.TRUE.,RIGIDBODY,DEBUG,RMAT)
!            ELSEIF (PERMDIST) THEN

!
! LOCAL ALIGNMENT TO IDENTIFY OPTIMAL "RIGID" GROUPS WITHIN THE GIVEN TOLERANCE.
!
             IF (LOCALPERMDIST) THEN
                CALL MINPERMRB(MI(I)%DATA%X,MI(J)%DATA%X,NATOMS, &
  &                              DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,D,DIST2,RIGIDBODY,RMAT)
             ELSE
                CALL MINPERMDIST(MI(I)%DATA%X,MI(J)%DATA%X,NATOMS, &
  &                              DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,D,DIST2,RIGIDBODY,RMAT)
             ENDIF
             IF (D < GEOMDIFFTOL) THEN
                  WRITE(*,'(3(A,G20.10))') ' CHECKPAIR> DISTANCE ',D,' IS LESS THAN TOLERANCE ',GEOMDIFFTOL, &
  &                                        ' - THIS SHOULD NOT HAPPEN'
                  CALL TSUMMARY
                  STOP
             ENDIF
          ELSE
!            IF (BULKT) THEN
!               CALL NEWMINDIST(MI(I)%DATA%X,MI(J)%DATA%X,NATOMS,D,BULKT,TWOD,"AX   ",.FALSE.,RIGIDBODY,DEBUG,RMAT)
!            ELSEIF (PERMDIST) THEN
!
! LOCAL ALIGNMENT TO IDENTIFY OPTIMAL "RIGID" GROUPS WITHIN THE GIVEN TOLERANCE.
!
             IF (LOCALPERMDIST) THEN
                CALL MINPERMRB(MI(I)%DATA%X,MI(J)%DATA%X,NATOMS, &
  &                              DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,D,DIST2,RIGIDBODY,RMAT)
             ELSE
                CALL MINPERMDIST(MI(I)%DATA%X,MI(J)%DATA%X,NATOMS, &
  &                              DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,D,DIST2,RIGIDBODY,RMAT)
             ENDIF
          ENDIF

     END SUBROUTINE CHECKPAIR

     FUNCTION GETDISTANCE(I,J)
          USE KEY,ONLY: DEBUG
          IMPLICIT NONE

          INTEGER,INTENT(IN) :: I,J
          DOUBLE PRECISION :: GETDISTANCE

          IF (I<J) THEN
               GETDISTANCE=MI(J)%DATA%D(I)
          ELSE IF (I>J) THEN
               GETDISTANCE=MI(I)%DATA%D(J)
          ELSE
               IF (DEBUG) PRINT *, 'GETDISTANCE> WARNING: I = J =',I
               GETDISTANCE=0.0D0
          ENDIF
     END FUNCTION GETDISTANCE

     FUNCTION GETINTERP(I,J)
          USE KEY,ONLY: DEBUG
          IMPLICIT NONE

          INTEGER,INTENT(IN) :: I,J
          DOUBLE PRECISION :: GETINTERP

          IF (I<J) THEN
               GETINTERP=MI(J)%DATA%INTERP(I)
          ELSE IF (I>J) THEN
               GETINTERP=MI(I)%DATA%INTERP(J)
          ELSE
               IF (DEBUG) PRINT *, 'GETINTERP> WARNING: I = J =',I
               GETINTERP=0.0D0
          ENDIF
     END FUNCTION GETINTERP

     SUBROUTINE SETDISTANCE(I,J,D)
          USE KEY,ONLY: DEBUG
          IMPLICIT NONE
          INTEGER,INTENT(IN) :: I,J
          DOUBLE PRECISION,INTENT(IN) :: D

          IF (I<J) THEN
               MI(J)%DATA%D(I)=D
          ELSE IF (I>J) THEN
               MI(I)%DATA%D(J)=D
          ELSE
               IF (DEBUG) PRINT *, 'SETDISTANCE> WARNING: I = J =',I
          ENDIF
     END SUBROUTINE SETDISTANCE

     SUBROUTINE SETINTERP(I,J,INTVALUE)
          USE KEY,ONLY: DEBUG
          IMPLICIT NONE
          INTEGER,INTENT(IN) :: I,J
          DOUBLE PRECISION,INTENT(IN) :: INTVALUE

          IF (I<J) THEN
               MI(J)%DATA%INTERP(I)=INTVALUE
          ELSE IF (I>J) THEN
               MI(I)%DATA%INTERP(J)=INTVALUE
          ELSE
               IF (DEBUG) PRINT *, 'SETINTERP> WARNING: I = J =',I
          ENDIF
     END SUBROUTINE SETINTERP

     DOUBLE PRECISION FUNCTION INTERPVALUE(Q1,Q2,DISTANCE) ! FINDS THE HIGHEST POINT ON A DISCRITISED LINEAR INTERPOLATION
          USE KEYDECIDE,ONLY : INTERPDIFF
          USE KEY,ONLY : DEBUG
          IMPLICIT NONE 
          DOUBLE PRECISION Q1(NOPT), Q2(NOPT), DELTAX(NOPT), DISTANCE, LOCALPOINTS(NOPT), MAXE, ENERGY, VNEW(NOPT)
          INTEGER NSTEPS, I

          DELTAX(1:NOPT) = ( Q1(1:NOPT) - Q2(1:NOPT) )*INTERPDIFF/DISTANCE
          NSTEPS=INT(DISTANCE/INTERPDIFF)
          MAXE=-1.0D100
          DO I=1,NSTEPS
             LOCALPOINTS(1:NOPT) = Q2(1:NOPT) + DELTAX*(I-1)
             CALL POTENTIAL(LOCALPOINTS,ENERGY,VNEW,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
             IF (ENERGY.GT.MAXE) MAXE=ENERGY
!            PRINT '(A,3G20.10)',' INTERPVALUE> I,ENERGY,MAXE=',I,ENERGY,MAXE
             IF (DEBUG) PRINT '(A,3G20.10)',' INTERPVALUE> I,ENERGY,MAXE=',I,ENERGY,MAXE
          ENDDO
          INTERPVALUE=MAXE

     END FUNCTION INTERPVALUE

END MODULE CONNECTUTILS
