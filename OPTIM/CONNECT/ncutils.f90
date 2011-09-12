!   CONNECT module is an implementation of a connection algorithm for finding rearrangement pathways.
!   Copyright (C) 2003-2006 Semen A. Trygubenko and David J. Wales
!   This file is part of CONNECT module. CONNECT module is part of OPTIM.
!
!   OPTIM is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   OPTIM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
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
!print *, "Energy of the TS found is the same as of the TS #",i,"; checking distance..."
                  CALL MINPERMDIST(TSTOCHECK%COORD,TS(I)%DATA%X, NATOMS, &
  &                                DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,D,DIST2,RIGIDBODY,RMAT)
!print *, "The distance is",d
                    IF (D<GEOMDIFFTOL) THEN
                         !print *, "This TS is already known."
                         ISNEWTS=.FALSE.
!                        IF (I > NTSOLD) THEN ! TWO OR MORE TS RETURNED FROM _ONE_ NEB RUN ARE THE SAME
!                             PRINT '(A,I8)', ' isnewts> Shortcut was found - TS is same as already saved TS #',i 
!                        ENDIF
                         PRINT '(A,I8,A,F20.10)', ' isnewts> transition state is same as number ',I,' energy=',TSTOCHECK%E
                         RETURN
                    ENDIF
               ENDIF
          ENDDO

          ISNEWTS=.TRUE.
     END FUNCTION ISNEWTS
    
     SUBROUTINE ISNEWMIN(E,COORD,MINPOS,NEW,REDOPATH,PERMUTE,INVERT,INDEX,I)
          USE KEYUTILS
          USE KEY,ONLY : DEBUG, RIGIDBODY, BULKT, TWOD, PERMDIST, AMBERT !msb50
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
!            IF (DEBUG) PRINT '(A,I6,2G20.10)',' isnewmin> I,E,MI(I)%DATA%E=',I,E,MI(I)%DATA%E
             IF (ABS(E-MI(I)%DATA%E) < EDIFFTOL) THEN
                CALL MINPERMDIST(COORD,MI(I)%DATA%X(1:3*NATOMS), NATOMS, &
  &                                DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,D,DIST2,RIGIDBODY,RMAT)
!               IF (DEBUG) PRINT '(A,G20.10)','isnewmin> minimum distance=',D
                IF (D<GEOMDIFFTOL) THEN
                   NEW=.FALSE.
                   MINPOS=I
!                  IF (DEBUG) PRINT '(A,I6)','isnewmin> MINPOS=',I
                   RETURN
                ENDIF
             ENDIF
          ENDDO

!                 SUCCESS=.FALSE. 
! !
! !  If they are the same minimum then GETPERM should be more reliable than mindist. DJW
! !  Even if the permutation-inversion is the identity!
! !  GETPERM changes the first argument, but not the second.
! !  However! GETPERM will fail if the stationary points are within EDiffTol but
! !  don;t line up sufficiently well. 
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
!                    ELSE ! don;t return unless PERMDIST is true - we may need to check the other minima for an exact match
!                       IF (INVERT.EQ.1) PRINT '(A,I6,A)',   &
!                              ' isnewmin> Permutational isomer of minimum ',i,' identified'
!                       IF (INVERT.EQ.-1) PRINT '(A,I6,A)',  &
!   &                          ' isnewmin> Permutation-inversion isomer of minimum ',i,' identified'
!                       IF (PERMDIST) THEN
!                          MINPOS=I   ! we should not need to save any permutation-inversion isomers if PERMDIST is true
!                          NEW=.FALSE.
!                          RETURN
!                       ENDIF
!                    ENDIF
!                 ENDIF
!              ENDIF
!           ENDDO
! 
!           IF ((MINPOS.GT.0).AND.REDOPATH) THEN
!              NEW=.FALSE. ! MINPOS was set in the last match
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
!                   print '(A)','calling mindist 3'
              CALL NEWMINDIST(COORD,COORD2,NATOMS,D,BULKT,TWOD,'AX   ',.True.,RIGIDBODY,DEBUG,RMAT)
              IF (D < GEOMDIFFTOL) THEN
                  SAME = .TRUE.
                  WRITE(*,'(a)') ' Same minimum found on either side of the pathway - rejecting this min--sad--min triple'
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
!         PRINT '(A,I6)',' addnewmin> NMIN=',NMIN
!         DO I=1,NMIN
!           PRINT '(A,I6,F20.10)',' addnewmin> I,energy=',I,MI(I)%DATA%E
!         ENDDO
          ALLOCATE( MI(NMIN)%DATA%D(NMIN-1),MI(NMIN)%DATA%NTRIES(NMIN-1),MI(NMIN)%DATA%C(1) )
          IF (INTERPCOSTFUNCTION) ALLOCATE( MI(NMIN)%DATA%INTERP(NMIN-1))

          DO I=1,NMIN-1
               CALL MINPERMDIST(MI(NMIN)%DATA%X(:), MI(I)%DATA%X(:), NATOMS, &
  &                             DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,D,DIST2,RIGIDBODY,RMAT)
!              IF (DEBUG) PRINT '(A,G20.10,2(A,I6))',' addnewmin> distance from MINPERMDIST/NEWMINDIST=',D,' for ',I,' and ',NMIN
               MI(NMIN)%DATA%D(I)=D  !  PASSING MIN(NMIN)%DATA%D(I) INSTEAD OF D DOES NOT WORK
                                     !  perhaps because min is an intrisic function!
               IF (INTERPCOSTFUNCTION) THEN
                  IF (INTCONSTRAINTT) THEN
                     MINCOORDS(1,1:NOPT)=MI(NMIN)%DATA%X(1:NOPT)
                     MINCOORDS(2,1:NOPT)=MI(I)%DATA%X(1:NOPT)
                     CALL MAKE_CONPOT(2,MINCOORDS)
!
! NMAXINT and NMININT are returned.
! 
                     FREEZENODESTLOCAL=FREEZENODEST
                     FREEZENODEST=.FALSE.
                     XYZLOCAL(1:NOPT)=MINCOORDS(1,1:NOPT)
                     XYZLOCAL(NOPT+1:2*NOPT)=MINCOORDS(2,1:NOPT)
                     INTIMAGESAVE=INTIMAGE
                     INTIMAGE=0
                     CALL CONGRAD2(NMAXINT,NMININT,CONSTRAINTE,XYZLOCAL,GGGLOCAL,EEELOCAL,IMGFREEZELOCAL,RMSLOCAL)
                     MI(NMIN)%DATA%INTERP(I)=CONSTRAINTE/2.0D0 ! energy per image
                     IF (DEBUG) PRINT '(A,I6,A,I6,2(A,G20.10))',' addnewmin> interpolation metric for ',NMIN,&
  &                  ' and ',I,' is ',MI(NMIN)%DATA%INTERP(I),' dist=',MI(NMIN)%DATA%D(I)
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
                     MI(NMIN)%DATA%INTERP(I)=CONSTRAINTE/2.0D0 ! energy per image
                     IF (DEBUG) PRINT '(A,I6,A,I6,2(A,G20.10))',' addnewmin> interpolation metric for ',NMIN,&
  &                  ' and ',I,' is ',MI(NMIN)%DATA%INTERP(I),' dist=',MI(NMIN)%DATA%D(I)
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
          OPEN(UNIT=11,FILE="ts.xyz")
          DO I=1,NTS !FOUND 
               WRITE(11,'(i3)') nopt/3 
               WRITE(11,'(a,i3)') 'ts #', i
               DO J=1,NOPT,3
                    WRITE(11,'(a5,1x,3f20.10)') ZSYM((j+2)/3),ts(i)%data%X(j),ts(i)%data%X(j+1),ts(i)%data%X(j+2)
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
               
               !connected ts
               P => MI(I)%DATA%CTS
               NULLIFY(MI(I)%DATA%CTS)
               ALLOCATE( MI(I)%DATA%CTS(EXTENT+1) )
               MI(I)%DATA%CTS(:EXTENT) = P
               MI(I)%DATA%CTS(EXTENT+1)= T
               DEALLOCATE(P)

               !connected minima
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
          
          ! storing new connection for ts and minima involved
          TS(K)%DATA%P=I
          TS(K)%DATA%M=J
          CALL ADDNEWCONNECTION(I,K,J)
          CALL ADDNEWCONNECTION(J,K,I)

          ! setting flags:
          IF (.NOT.MI(I)%DATA%S .AND. .NOT.MI(I)%DATA%F .AND. .NOT.MI(J)%DATA%S .AND. .NOT.MI(J)%DATA%F) THEN
               PRINT *, 'Connection established between members of the U set.'
               RETURN 
          ELSEIF ( (MI(I)%DATA%S.AND.MI(J)%DATA%S).OR.(MI(I)%DATA%F.AND.MI(J)%DATA%F) ) THEN
               IF (MI(I)%DATA%S) THEN
                    PRINT *, 'Alternative path found between members of the S set.'
               ELSE
                    PRINT *, 'Alternative path found between members of the F set.'
               ENDIF
               RETURN
          ELSEIF ( (MI(I)%DATA%S.AND.MI(J)%DATA%F).OR.(MI(J)%DATA%S.AND.MI(I)%DATA%F) ) THEN
               IF (FINISHED) RETURN ! Must not call GETPATHWAY twice - segmentation fault! DJW
               FINISHED     =.TRUE.
               CALL GETPATHWAY(I,J)
               RETURN
          ELSEIF ( (.NOT.MI(I)%DATA%S .AND. .NOT.MI(I)%DATA%F .AND. MI(J)%DATA%S) .OR. &
                 & (.NOT.MI(J)%DATA%S .AND. .NOT.MI(J)%DATA%F .AND. MI(I)%DATA%S) ) THEN
               IF (MI(J)%DATA%S) THEN
                    WRITE(CHR,'(i7)') i
                    CALL UPDATELINK(I,J)
               ELSE
                    WRITE(CHR,'(i7)') j
                    CALL UPDATELINK(J,I)
               ENDIF
               WRITE(*,'(1x,a)') 'Unconnected minimum '//trim(adjustl(chr))//' found its way to S set.'
          ELSEIF ( (.NOT.MI(I)%DATA%S .AND. .NOT.MI(I)%DATA%F .AND. MI(J)%DATA%F) .OR. &
                 & (.NOT.MI(J)%DATA%S .AND. .NOT.MI(J)%DATA%F .AND. MI(I)%DATA%F) ) THEN
               IF (MI(J)%DATA%F) THEN
                    WRITE(CHR,'(i7)') i
                    CALL UPDATELINK(I,J)
               ELSE
                    WRITE(CHR,'(i7)') j
                    CALL UPDATELINK(J,I)
               ENDIF
               WRITE(*,'(1x,a)') 'Unconnected minimum '//trim(adjustl(chr))//' found its way to F set.'
          ENDIF 
     END SUBROUTINE NEWCONNECTION

     ! Semen Срд Лип 19 14:37:58 BST 2006: updatelink is a new version of updatelinkold;
     ! should be faster, easier to maintain and debug; updatelinkold uses the same data structures
     ! so should still work --- to enable swap updatelink and updatelinkold.
     RECURSIVE SUBROUTINE UPDATELINK(I,J)
          use ConnectData,only:listlength
          IMPLICIT NONE
          INTEGER,INTENT(IN) :: I,J ! minima i should be brought from u set to set to which belongs minimum j;
          INTEGER :: K,M

          DEPTH = DEPTH+1 ! depth level of the recursion
!         PRINT '(A,3I8)','UPDATELINK> I,J,DEPTH=',I,J,DEPTH
       
          ! add j to the list of minima to be avoided
          IF (listlength==0) THEN 
               ALLOCATE(START2)
               NULLIFY(START2%NEXT)
               DUMMY2 => START2
          ELSE
               ALLOCATE(DUMMY2%NEXT)
               DUMMY2 => DUMMY2%NEXT
               NULLIFY(DUMMY2%NEXT)
          ENDIF
          listlength=listlength+1
          DUMMY2%I=J

          ! update flags
          MI(I)%DATA%S = MI(J)%DATA%S
          MI(I)%DATA%F = MI(J)%DATA%F
          !write(*,'(a)',advance='no') '.'

          ! update all connections of minimum i
       U: DO K=1,MI(I)%DATA%C(1)

               ! do not try connections to minima in the list
               DUMMY2 => START2
               IF (MI(I)%DATA%CMIN(K)==DUMMY2%I) CYCLE U
               DO M=2,listlength
                    DUMMY2 => DUMMY2%NEXT
                    IF (MI(I)%DATA%CMIN(K)==DUMMY2%I) CYCLE U
               ENDDO

               CALL UPDATELINK(MI(I)%DATA%CMIN(K),I)
          ENDDO U
          
          if (depth==1) then
               do m=1,listlength
                    dummy2=>start2%next
                    nullify(start2%next)
                    deallocate(start2)
                    start2=>dummy2
               enddo
               listlength=0
          endif

          depth = depth-1
     END SUBROUTINE UPDATELINK

     RECURSIVE SUBROUTINE UPDATELINKOLD(I,J)
          IMPLICIT NONE
          INTEGER,INTENT(IN)          :: I,J ! MINIMA I SHOULD BE BROUGHT FROM U SET TO SET TO WHICH BELONGS MINIMUM J;
                                             ! all connections must be updated; connection to minimum j is ignored
          INTEGER :: K,M
       
          DEPTH = DEPTH+1 ! DEPTH LEVEL OF THE RECURSION; EQUALS TO THE LENGTH OF THE LIST OF MINIMA TO BE IGNORED
!         PRINT '(A,3I8)','UPDATELINK> I,J,DEPTH=',I,J,DEPTH

          ! add j to the list of minima to be avoided
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

          ! update flags
          MI(I)%DATA%S = MI(J)%DATA%S
          MI(I)%DATA%F = MI(J)%DATA%F
          !write(*,'(a)',advance='no') '.'

          ! update all connections of minimum i
       U: DO K=1,MI(I)%DATA%C(1)

               ! do not try connections to minima in the list
               DUMMY2 => START2
               IF (MI(I)%DATA%CMIN(K)==DUMMY2%I) CYCLE
               DO M=2,DEPTH
                    DUMMY2 => DUMMY2%NEXT
                    IF (MI(I)%DATA%CMIN(K)==DUMMY2%I) CYCLE U
               ENDDO

               CALL UPDATELINK(MI(I)%DATA%CMIN(K),I)
          ENDDO U
          
          ! remove last minimum from the list
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
                                    ! yet i belongs to S, j - to F (or vice versa)
          ! constructing the chain
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
! Without this line ifort debug complains when we hit a TSINDEX that is not assigned
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
! Note for STOCKT. The path.<n>.xyz files have got the centre of mass coordinates plus a 
! unit vector. We need to rotate all of these appropriately as if each one represents
! an atom, i.e. using NEWROTGEOM, not NEWROTGEOMSTOCK.
!
! Remaining problems for STOCK:
! The energy is invariant to a change in sign of all the dipoles.
! Aligning paths properly might therefore need an inversion of all
! the dipole moments. We are actually aligning the path.<n>.xyz files,
! where the dipoles are present as extra unit vectors which are treated
! as atoms in the following subroutine. 
! A good solution would probably be to define an unambiguous sense for the
! dipoles so that this problem doesn't arise. It could also cause problems
! in recognising identical minima.
!
! Merges path output files to produce full pathway for the rearrangement;
! frames in path are reversed as needed;
!
     SUBROUTINE MERGEXYZEOFS  
          USE KEY, ONLY: FILTH,UNRST,FILTHSTR,DEBUG,RIGIDBODY,BULKT,TWOD,STOCKT,RBAAT,PERMDIST, AMHT,SEQ,NTSITES
          USE KEYUTILS        ! frames in bits that are glued together are rotated accordingly;
          USE COMMONS,ONLY : PARAM1,PARAM2,PARAM3
          USE AMHGLOBALS, ONLY : NMRES

          IMPLICIT NONE       ! prerequisites: chain of min/ts constructed; assumes path is dumping plus side of the path
                              ! first, and there are no blank lines after last frame (!)
                              ! does somewhat similar with EofS.ts files as well..
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
             IF (RBAAT) OPEN(UNIT=55,FILE='rbpath.xyz',status='UNKNOWN')
!            OPEN(UNIT=50,FILE='path.xyz',status='replace',action='write')
             OPEN(UNIT=50,FILE='path.xyz',status='UNKNOWN')
             IF (UNRST) OPEN(UNIT=51,FILE='path.unr.xyz',status='replace',action='write')
!            OPEN(UNIT=41,FILE='EofS',status='replace',action='write')
             OPEN(UNIT=41,FILE='EofS',status='UNKNOWN')
          ELSE
             IF (RBAAT) THEN
                WRITE(RBPSTRING,'(A)') 'rbpath.xyz.'//TRIM(ADJUSTL(FILTHSTR))
!               OPEN(UNIT=55,FILE=PSTRING,STATUS='replace',action='write')
                OPEN(UNIT=55,FILE=PSTRING,STATUS='UNKNOWN')
             ENDIF
             WRITE(PSTRING,'(A)') 'path.xyz.'//TRIM(ADJUSTL(FILTHSTR))
             IF (UNRST) WRITE(PSTRING2,'(A)') 'path.unr.xyz.'//TRIM(ADJUSTL(FILTHSTR))
             WRITE(ESTRING,'(A)') 'EofS.'//TRIM(ADJUSTL(FILTHSTR))
!            OPEN(UNIT=50,FILE=PSTRING,STATUS='replace',action='write')
             OPEN(UNIT=50,FILE=PSTRING,STATUS='UNKNOWN')
             IF (UNRST) OPEN(UNIT=51,FILE=PSTRING2,STATUS='replace',action='write')
!            OPEN(UNIT=41,FILE=ESTRING,STATUS='replace',action='write')
             OPEN(UNIT=41,FILE=ESTRING,STATUS='UNKNOWN')
          ENDIF
          
          DUMMY => START
          SLAST=0.0D0
          I=1      ! DC430 > I loops over paths
          DO
             IF (.NOT.ASSOCIATED(DUMMY%NEXT)) EXIT
             REVERSEORDER=.TRUE.
             IF (TS( MI(DUMMY%I)%DATA%CTS(DUMMY%J) )%DATA%P == DUMMY%I) REVERSEORDER=.FALSE.
             CALL MKFNAMES(MI(DUMMY%I)%DATA%CTS(DUMMY%J),FILTH,FILTHSTR,ITSTRING,EOFSSTRING)
!
!  If READSP is true then the pathway coordinates and energy files may not be present.
!  We could either ignore this and skip on, or attempt to recalculate the pathway.
!  If we recalculate the pathway, there is a danger that it might not link the expected
!  minima!
!
             INQUIRE(FILE=ITSTRING,EXIST=EXTEST)
             IF (.NOT.EXTEST) THEN
                REVERSEORDER=.TRUE.
                INDEXTS=MI(DUMMY%I)%DATA%CTS(DUMMY%J)
                PRINT '(A,I6,A)',' pathway files for ts ',INDEXTS,' are missing, attempting to recreate them'
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
                      PRINT '(A)',' energy from plus side of path matches expected plus minimum'

                      CALL MINPERMDIST(MI(TS(INDEXTS)%DATA%P)%DATA%X,QPLUS,NATOMS,DEBUG, &
    &                         PARAM1,PARAM2,PARAM3,BULKT,TWOD,DIST,DIST2,RIGIDBODY,RMAT)
                      IF (D.LE.GEOMDIFFTOL) THEN
                         PRINT '(A)',' geometry from plus side of path matches expected plus minimum'
                         IF (ABS(EMINUS-MI(TS(INDEXTS)%DATA%M)%DATA%E).LE.EDIFFTOL) THEN
                            PRINT '(A)',' energy from minus side of path matches expected minus minimum'
                            CALL MINPERMDIST(MI(TS(INDEXTS)%DATA%M)%DATA%X,QMINUS,NATOMS,DEBUG, &
    &                                           PARAM1,PARAM2,PARAM3,BULKT,TWOD,DIST,DIST2,RIGIDBODY,RMAT)
                            IF (D.LE.GEOMDIFFTOL) THEN
                               PRINT '(A)',' geometry from minus side of path matches expected minus minimum'
                               IF (TS(INDEXTS)%DATA%P == DUMMY%I) REVERSEORDER=.FALSE.
                               PATHFAILT=.FALSE.
                            ELSE
                               PRINT '(A)',' geometry from minus side of path does not match expected minus minimum'
                            ENDIF
                         ELSE
                            PRINT '(A)',' energy from minus side of path does not match expected minus minimum'
                         ENDIF
                      ELSE
                         PRINT '(A)',' geometry from plus side of path does not match that expected for the plus minimum'
                      ENDIF
                   ENDIF
                   IF (PATHFAILT) THEN
                      IF (ABS(EPLUS-MI(TS(INDEXTS)%DATA%M)%DATA%E).LE.EDIFFTOL) THEN
                         PRINT '(A)',' energy from plus side of path matches expected minus minimum'
                         CALL MINPERMDIST(MI(TS(INDEXTS)%DATA%M)%DATA%X,QPLUS,NATOMS,DEBUG, &
    &                                        PARAM1,PARAM2,PARAM3,BULKT,TWOD,DIST,DIST2,RIGIDBODY,RMAT)
                         IF (D.LE.GEOMDIFFTOL) THEN
                            PRINT '(A)',' geometry from plus side of path matches expected minus minimum'
                            IF (ABS(EMINUS-MI(TS(INDEXTS)%DATA%P)%DATA%E).LE.EDIFFTOL) THEN
                               PRINT '(A)',' energy from minus side of path matches expected plus minimum'
                               CALL MINPERMDIST(MI(TS(INDEXTS)%DATA%P)%DATA%X,QMINUS,NATOMS,DEBUG, &
    &                                              PARAM1,PARAM2,PARAM3,BULKT,TWOD,DIST,DIST2,RIGIDBODY,RMAT)
                               IF (D.LE.GEOMDIFFTOL) THEN
                                  PRINT '(A)',' geometry from minus side of path matches expected plus minimum'
                                  PATHFAILT=.FALSE.
                                  IF (TS(INDEXTS)%DATA%M == DUMMY%I) REVERSEORDER=.FALSE.
                               ELSE
                                  PRINT '(A)',' geometry from minus side of path does not match expected plus minimum'
                               ENDIF
                            ELSE
                               PRINT '(A)',' energy from minus side of path does not match expected plus minimum'
                            ENDIF
                         ELSE
                            PRINT '(A)',' geometry from plus side of path does not match that expected for the minus minimum'
                         ENDIF
                      ENDIF
                   ENDIF
                ENDIF
                IF (PATHFAILT) PRINT '(A)',' ERROR - attempt to fill in pathway is not consistent with saved data'
             ENDIF
               
             !XYZ BIT

             !reading part of the path that corresponds to ts dummy%j
             OPEN(UNIT=40,FILE=ITSTRING,STATUS='old',action="read")
             K=1   ! DC430 > K loops over frames in a given path corresponding to a single TS
             DO
                  READ(40,*,IOSTAT=EOF) !NEVIDOMO SKIL'ky fraimiv
                  IF (EOF==0) THEN
                       BACKSPACE 40
                  ELSE
                       EXIT
                  ENDIF
                  IF (K<0) PRINT*,'k=',k ! stupid fix for stupid Sun compiler bug
                                         ! need to access k to prevent SEGV       DJW
                    
                  ! stvoryty novyi element
                  ALLOCATE(TAIL); NULLIFY(TAIL%NEXT,TAIL%PREVIOUS,TAIL%Q,TAIL%SYM,TAIL%LINE)
                    
                  ! i zapovnyty yogo danymy
                  ALLOCATE(TAIL%Q(3*NATOMS),TAIL%SYM(NATOMS))
                  READ(40,'(a)')
                  READ(40,'(a)') tail%comment
                  IF (STOCKT) THEN
                     DO J=1,(NATOMS/2)
                        READ(40,'(a5,1x,3f20.10,tr13,3f20.10)') tail%sym(j), &
  &                     tail%q(3*(j-1)+1),tail%q(3*(j-1)+2),tail%q(3*(j-1)+3), &
  &                     tail%q(3*((NATOMS/2)+j-1)+1),tail%q(3*((NATOMS/2)+j-1)+2),tail%q(3*((NATOMS/2)+j-1)+3)
                        DLENGTH=TAIL%Q(3*((NATOMS/2)+J-1)+1)**2 &
  &                            +TAIL%Q(3*((NATOMS/2)+J-1)+2)**2 & 
  &                            +TAIL%Q(3*((NATOMS/2)+J-1)+3)**2 
                        IF (ABS(DLENGTH-1.0D0).GT.1.0D-5) THEN
                           PRINT '(A,I8,G20.10)','1 merge> J,length=',J,DLENGTH
                           PRINT '(A,3G20.10)','unit vector: ',TAIL%Q(3*((NATOMS/2)+J-1)+1:3*((NATOMS/2)+J-1)+3)
                           STOP
                        ENDIF
                     ENDDO
                  ELSE IF(RBAAT) THEN
                     DO J=1,(NATOMS/2)
                        READ(40,'(a5,1x,6f20.10)') tail%sym(j), &
  &                     tail%q(3*(j-1)+1),tail%q(3*(j-1)+2),tail%q(3*(j-1)+3), &
  &                     tail%q(3*((NATOMS/2)+j-1)+1),tail%q(3*((NATOMS/2)+j-1)+2),tail%q(3*((NATOMS/2)+j-1)+3)
                     ENDDO
                  ELSEIF (AMHT) THEN
                     GLY_COUNT = 0
                     NDUMMY=0
                     DO J=1,NMRES
                        IF (SEQ(J).EQ.8) THEN
                           READ(40,'(a5,1x,3f20.10)')AMHDUMMY,XDUMMY,XDUMMY,XDUMMY
                           NDUMMY=NDUMMY+1
                           READ(40,'(a5,1x,3f20.10)') tail%sym(NDUMMY),tail%Q(3*(NDUMMY-1)+1),tail%Q(3*(NDUMMY-1)+2),  &
     &                                                tail%Q(3*(NDUMMY-1)+3)
                           NDUMMY=NDUMMY+1
                           READ(40,'(a5,1x,3f20.10)') tail%sym(NDUMMY),tail%Q(3*(NDUMMY-1)+1),tail%Q(3*(NDUMMY-1)+2),  &
     &                                                tail%Q(3*(NDUMMY-1)+3)
                        ELSE
                           NDUMMY=NDUMMY+1
                           READ(40,'(a5,1x,3f20.10)') tail%sym(NDUMMY),tail%Q(3*(NDUMMY-1)+1),tail%Q(3*(NDUMMY-1)+2),  &
     &                                                tail%Q(3*(NDUMMY-1)+3)
                           NDUMMY=NDUMMY+1
                           READ(40,'(a5,1x,3f20.10)') tail%sym(NDUMMY),tail%Q(3*(NDUMMY-1)+1),tail%Q(3*(NDUMMY-1)+2),  &
     &                                                tail%Q(3*(NDUMMY-1)+3)
                           NDUMMY=NDUMMY+1
                           READ(40,'(a5,1x,3f20.10)') tail%sym(NDUMMY),tail%Q(3*(NDUMMY-1)+1),tail%Q(3*(NDUMMY-1)+2),  &
     &                                                tail%Q(3*(NDUMMY-1)+3)
                        ENDIF
                     ENDDO
                  ELSE
                     DO J=1,NATOMS
                        READ(40,'(a5,1x,3f20.10)') tail%sym(j),tail%q(3*(j-1)+1),tail%q(3*(j-1)+2),tail%q(3*(j-1)+3)
                     ENDDO
                  ENDIF
                  
                  IF (ASSOCIATED(TMP)) THEN ! ODYN ABO DEKIL'ka lanok vge isnye; teper tmp mistyt' poperednu lanku
                       TMP%NEXT => TAIL
                       TAIL%PREVIOUS => TMP ! ZVOROTNIY ZVYAZOK
                  ENDIF

                  TMP=>TAIL 
                  IF (K==1) HEAD => TAIL
                  K=K+1
             ENDDO
             CLOSE(40)
               
!  Writing frames to path.xyz file in correct order + rotated as needed
!  For Stockmayer the coordinates read in from path.*.xyz contain the
!  vector components, rather than than spherical polar angles. Hence
!  we want to call NEWROTGEOM to rotate all of the sites.
!  We also need STOCKT turned off in MINPERMDIST, which is why we
!  use a LOCALSTOCK in this subroutine.

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
!  May need to invert dipoles for STOCK.
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
                  !open(unit=40,file=EofSString,status='old',action='read',position='append')
                  ! PGI complains about read and append together ... hm...
                  OPEN(UNIT=40,FILE=EOFSSTRING,STATUS='old',position='append')
                  BACKSPACE 40
                  READ(40,'(f20.10)') Sfirst
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
                       READ(40,'(a)') head%line(1)
                       TMP=>HEAD
                  ENDDO
                  DO
                       READ(HEAD%LINE(1),'(2f20.10)') s,e
                       DEALLOCATE(HEAD%LINE)
                       WRITE(UNIT=41,FMT='(2f20.10)') Sfirst - s,e
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
                  OPEN(UNIT=40,FILE=EOFSSTRING,STATUS='old',action='read')
                  J=1
                  DO
                     READ(UNIT=40,FMT='(2f20.10)',iostat=eof) s,e
                     IF (.NOT.EOF==0) EXIT
                       IF (J==1) SFIRST=ABS(S) + SLAST
                       WRITE(UNIT=41,FMT='(2f20.10)') Sfirst + s,e
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
             WRITE(50,'(i6)') (Natoms/2)
             WRITE(50,'(1x,a)') trim(adjustl(c))
             DO J=1,(NATOMS/2)
!               WRITE(50,'(a5,1x,3f20.10,a13,3f20.10)') s(j),q(3*(j-1)+1),q(3*(j-1)+2),q(3*(j-1)+3), &
! &                             ' atom_vector ', &
! &                             q(3*((NATOMS/2)+j-1)+1),q(3*((NATOMS/2)+j-1)+2),q(3*((NATOMS/2)+j-1)+3)
                WRITE(50,'(A5,1X,6F20.10)') S(j),Q(3*(j-1)+1),Q(3*(j-1)+2),Q(3*(j-1)+3), &
  &                             q(3*((NATOMS/2)+j-1)+1),q(3*((NATOMS/2)+j-1)+2),q(3*((NATOMS/2)+j-1)+3)
             ENDDO
          ELSE IF (RBAAT) THEN
             WRITE(50,'(i6)') (Natoms/2)
             WRITE(50,'(1x,a)') trim(adjustl(c))
             WRITE(55,'(i6)') NTSITES
             WRITE(55,'(1x,a)') trim(adjustl(c))
             DO J=1,(NATOMS/2)
                WRITE(50,'(a5,1x,6f20.10)') s(j),q(3*(j-1)+1),q(3*(j-1)+2),q(3*(j-1)+3), &
  &                             q(3*((NATOMS/2)+j-1)+1),q(3*((NATOMS/2)+j-1)+2),q(3*((NATOMS/2)+j-1)+3)
             ENDDO
             CALL SITEPOS(Q,SITES)
             DO J=1,NTSITES
                WRITE(55,'(A5,1X,3F20.10)') 'LA ',SITES(3*(J-1)+1),SITES(3*(J-1)+2),SITES(3*(J-1)+3)
             ENDDO
          ELSE
             WRITE(50,'(i6)') Natoms
             WRITE(50,'(1x,a)') trim(adjustl(c))
             DO J=1,NATOMS
                WRITE(50,'(a5,1x,3f20.10)') s(j),q(3*(j-1)+1),q(3*(j-1)+2),q(3*(j-1)+3)
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

          WRITE(51,'(i6)') 2*Natoms-2 ! jmc for carbons plus dummy peptide atoms, want 2*Natoms-2...
          WRITE(51,'(1x,a)') trim(adjustl(c))
          DO J=1,NATOMS
               WRITE(51,'(a5,1x,3f20.10)') s(j),q(3*(j-1)+1),q(3*(j-1)+2),q(3*(j-1)+3)
          ENDDO
          
          DO K1=1,(NATOMS/2)-1
             DO K2=1,3
                PEPCOORDS(6*(K1-1)+K2)=(2.0D0*Q(6*(K1-1)+K2)+Q(6*K1+K2))/3.0D0
                PEPCOORDS(6*(K1-1)+K2+3)=(Q(6*(K1-1)+K2)+2.0D0*Q(6*K1+K2))/3.0D0
             ENDDO
          ENDDO
          WRITE(51,'(A2,4X,3F20.10)') ('O ',pepcoords(6*(K1-1)+1),pepcoords(6*(K1-1)+2),pepcoords(6*(K1-1)+3) &
                                       & ,K1=1,(NATOMS/2)-1)
          WRITE(51,'(A2,4X,3F20.10)') ('N ',pepcoords(6*(K1-1)+4),pepcoords(6*(K1-1)+5),pepcoords(6*(K1-1)+6) &
                                       & ,K1=1,(NATOMS/2)-1)

     END SUBROUTINE WRITEFRAMEUNRES

     SUBROUTINE REALLOCATETSRACK  ! RESIZES TS RACK: NEWSIZE = OLDSIZE*INCR
          IMPLICIT NONE
          
          INTEGER, PARAMETER :: INCR=2
          INTEGER :: I
          
          TEMPTSRACK => TS
          TSRACKSIZE=TSRACKSIZE*INCR
          WRITE(CHR,'(i5)') tsracksize
          PRINT '(1x,a)', 'TS rack redimentioned to hold '//trim(adjustl(chr))//' ts.'
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
          WRITE(CHR,'(i5)') minracksize
          PRINT '(1x,a)', 'Minima rack redimentioned to hold '//trim(adjustl(chr))//' minima.'
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
                    MYSTR='points'//trim(adjustl(mystr))//'.out'
                    OPEN(UNIT=38,FILE=TRIM(ADJUSTL(MYSTR)),STATUS='unknown',form='unformatted',access='direct',recl=reclen)
                    WRITE(38,REC=1) (MI(DUMMY%I)%DATA%X(I),I=1,3*NATOMS)
                    CLOSE(38)
                    IF (ASSOCIATED(DUMMY%NEXT)) THEN
                         NSP=NSP+1
                         WRITE(MYSTR,*) NSP
                         MYSTR='points'//trim(adjustl(mystr))//'.out'
                         OPEN(UNIT=38,FILE=TRIM(ADJUSTL(MYSTR)),STATUS='unknown',form='unformatted',access='direct',recl=reclen)
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
               OPEN(UNIT=38,FILE="POINTS.TS",STATUS='unknown',form='unformatted',access='direct',recl=reclen)
               OPEN(UNIT=39,FILE="TS.DATA",STATUS='unknown',form='formatted')
               OPEN(UNIT=40,FILE="POINTS.MIN",STATUS='unknown',form='unformatted',access='direct',recl=reclen)
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
                    WRITE(39,'(2i10)') rec1,rec2
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
!  File path.info is now opened in keyword.f
!  
!    IF (FILTH.EQ.0) THEN
!       WRITE(PINFOSTRING,'(A9)') 'path.info'
!    ELSE
!       WRITE(PINFOSTRING,'(A)') 'path.info.'//TRIM(ADJUSTL(FILTHSTR))
!    ENDIF
!    IF (MACHINE) THEN
!         OPEN(UNIT=88,FILE=PINFOSTRING,STATUS='UNKNOWN',form='unformatted')
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
!               WRITE(88,'(3F20.10)') mi(dummy%i)%data%E, DIHE, ALLANG
           WRITE(88,'(3F20.10)') mi(dummy%i)%data%E, DIHE
        ELSE
           WRITE(88,'(F20.10)') mi(dummy%i)%data%E
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
!CALL INT_FROM_CART(.TRUE.,.FALSE.)
!CALL CHAINBUILD
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
! jmc           WRITE(88,'(3G20.10)') (FSAVEMIN(J2,J1),J2=1,3*NATOMS)
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
!CALL INT_FROM_CART(.TRUE.,.FALSE.)
!CALL CHAINBUILD
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
                    IF (DIAG(J2).LT.0.0D0) PRINT *,'Higher order saddle found in pathway - ts ',i,'eigenvalue ',DIAG(J2)
                 END DO
                 WRITE(88,'(3G20.10)') (DIAG(J2),J2=1,3*NATOMS)
              ENDIF
           ELSEIF (AMHT) THEN
              PRINT '(A)','in ts AMH block writing 1 C1'
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
                 IF (DIAG(J2).LT.0.0D0) PRINT *,'Higher order saddle found in pathway - ts ',i,'eigenvalue ',DIAG(J2)
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
!  Dump the latest min-sad-min triple to path.info in the usual format
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
!  File path.info is now opened in keyword.f
!  
!    IF (FILTH.EQ.0) THEN
!       WRITE(PINFOSTRING,'(A9)') 'path.info'
!    ELSE
!       WRITE(PINFOSTRING,'(A)') 'path.info.'//TRIM(ADJUSTL(FILTHSTR))
!    ENDIF
!    if (machine) then
!         OPEN(UNIT=88,FILE=PINFOSTRING,STATUS='UNKNOWN',form='unformatted',POSITION='APPEND')
!    else
!         OPEN(UNIT=88,FILE=PINFOSTRING,STATUS='UNKNOWN',POSITION='APPEND')
!    endif

     TSFRQDONE=.FALSE.  ! ASSUME THAT WE NEVER KNOW THE FREQUENCIES
     MINFRQDONE=.FALSE. ! ASSUME THAT WE NEVER KNOW THE FREQUENCIES
!
!  First dump the + minimum.
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
!CALL INT_FROM_CART(.TRUE.,.FALSE.)
!CALL CHAINBUILD
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
! jmc   WRITE(88,'(3G20.10)') (FSAVEMIN(J2,J1),J2=1,3*NATOMS)
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
! now the transition state
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
!CALL INT_FROM_CART(.TRUE.,.FALSE.)
!CALL CHAINBUILD
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
!              OPEN(unit=955, file='hessdump', status='unknown')
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
                 IF (DIAG(J2).LT.0.0D0) PRINT *,'Higher order saddle found in pathway - ts ',i,'eigenvalue ',DIAG(J2)
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
!  Finally dump the - minimum.
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
!CALL INT_FROM_CART(.TRUE.,.FALSE.)
!CALL CHAINBUILD
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
! jmc   WRITE(88,'(3G20.10)') (FSAVEMIN(J2,J1),J2=1,3*NATOMS)
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

     KNOWH = .FALSE. ! needed otherwise the next TS search will use the wrong Hessian, if one is required.
     CALL FLUSH(88,ISTAT)

     END SUBROUTINE MAKEALLPATHINFO

     SUBROUTINE CHECKPAIR(I,J,PERMTEST) ! checks that minima I and J are different and puts them into the orientation with minimal D
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
             WRITE(*,'(/1x,a,2i5,a)') "Energies of the minima in the pair ",i,j," are the same - checking distance ..."

!            IF (BULKT) THEN
!               CALL NEWMINDIST(MI(I)%DATA%X,MI(J)%DATA%X,NATOMS,D,BULKT,TWOD,'AX   ',.True.,RIGIDBODY,DEBUG,RMAT)
!            ELSEIF (PERMDIST) THEN

!
! Local alignment to identify optimal "rigid" groups within the given tolerance.
!
             IF (LOCALPERMDIST) THEN
                CALL MINPERMRB(MI(I)%DATA%X,MI(J)%DATA%X,NATOMS, &
  &                              DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,D,DIST2,RIGIDBODY,RMAT)
             ELSE
                CALL MINPERMDIST(MI(I)%DATA%X,MI(J)%DATA%X,NATOMS, &
  &                              DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,D,DIST2,RIGIDBODY,RMAT)
             ENDIF
             IF (D < GEOMDIFFTOL) THEN
                  WRITE(*,'(3(A,G20.10))') ' checkpair> Distance ',D,' is less than tolerance ',GEOMDIFFTOL, &
  &                                        ' - this should not happen'
                  CALL TSUMMARY
                  STOP
             ENDIF
          ELSE
!            IF (BULKT) THEN
!               CALL NEWMINDIST(MI(I)%DATA%X,MI(J)%DATA%X,NATOMS,D,BULKT,TWOD,"AX   ",.False.,RIGIDBODY,DEBUG,RMAT)
!            ELSEIF (PERMDIST) THEN
!
! Local alignment to identify optimal "rigid" groups within the given tolerance.
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
               IF (DEBUG) PRINT *, 'GetDistance> WARNING: i = j =',i
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
               IF (DEBUG) PRINT *, 'getinterp> WARNING: i = j =',i
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
               IF (DEBUG) PRINT *, 'SetDistance> WARNING: i = j =',i
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
               IF (DEBUG) PRINT *, 'SetINTERP> WARNING: i = j =',i
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
!            PRINT '(A,3G20.10)',' interpvalue> i,energy,maxe=',i,energy,maxe
             IF (DEBUG) PRINT '(A,3G20.10)',' interpvalue> i,energy,maxe=',i,energy,maxe
          ENDDO
          INTERPVALUE=MAXE

     END FUNCTION INTERPVALUE

END MODULE CONNECTUTILS
