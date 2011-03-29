!
!     FOLLOWS THE SAME PROCEDURE AS IN MINPERMDIST.F90 TO ALIGN A SYSTEM OF RIGID BODIES BASED
!     ON THE CENTRES OF MASS ONLY.
!     ----------------------------------------------------------------------------------------------
 
      SUBROUTINE MINPERMDISTRBCOM(COORDSB,COORDSA,DISTANCE,DIST2,QBEST,RMATBEST,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT)
!     DISTANCE RETURNS THE SQUARED DISTANCE

      USE COMMONS, ONLY : NATOMS
      USE KEY,ONLY : NPERMGROUP, NPERMSIZE, PERMGROUP, NSETS, SETS, GEOMDIFFTOL, EFIELDT, &
     &               NRBGROUP

      IMPLICIT NONE

      INTEGER, PARAMETER :: MAXIMUMTRIES = 100
      INTEGER            :: NPERM, PATOMS, NTRIES, NSIZE, JMAX, LOCMAX(1), J1, J2, J3, INFO
      INTEGER            :: INVERT, NORBIT1, NORBIT2, PERM(NATOMS), NCHOOSE2, NDUMMY, LPERM(NATOMS), NCHOOSE1
      INTEGER            :: NEWPERM(NATOMS), ALLPERM(NATOMS), SAVEPERM(NATOMS)
      DOUBLE PRECISION   :: COORDSA(3*NATOMS), COORDSB(3*NATOMS), DISTANCE, DISTWP, DIST2, TEMPA(9*NATOMS) 
      DOUBLE PRECISION   :: DUMMYA(3*NATOMS), DUMMYB(3*NATOMS), DUMMY(3*NATOMS), DUMMYWP(3*NATOMS)
!      DOUBLE PRECISION   :: XA(3*NATOMS*NRBSITES/2),  XB(3*NATOMS*NRBSITES/2), XBS(3*NATOMS*NRBSITES/2)
      DOUBLE PRECISION   :: XA(3*NATOMS),  XB(3*NATOMS), XBS(3*NATOMS)
      DOUBLE PRECISION   :: XTMP(3*NATOMS)
      DOUBLE PRECISION   :: RMAT(3,3), RMATI(3,3), ENERGY, VNEW(3*NATOMS), DX, DY, DZ, RMS, DBEST, XBEST(3*NATOMS)
      DOUBLE PRECISION   :: ROTA(3,3), ROTINVA(3,3), ROTB(3,3), ROTBINV(3,3), RMATCUMUL(3,3), LMAT(3,3)
      DOUBLE PRECISION   :: ROTINVBBEST(3,3), ROTABEST(3,3), RMATBEST(3,3), RMATWP(3,3)
      DOUBLE PRECISION   :: CMAX, CMAY, CMAZ, CMBX, CMBY, CMBZ
      DOUBLE PRECISION   :: PDUMMYA(3*NATOMS), PDUMMYB(3*NATOMS), LDISTANCE, XDUMMY, BOXLX, BOXLY, BOXLZ, WORSTRAD
      DOUBLE PRECISION   :: Q(4), Q1(4), Q2(4), AMAT(4,4), BMAT(4,4), DIAG(4), P(3)
      DOUBLE PRECISION   :: ST, THETA, THETAH, FCT, DUMMYC(3*NATOMS), DUMMYD(3*NATOMS)
      LOGICAL            :: DEBUG, BULKT
      DOUBLE PRECISION   :: BESTA(3*NATOMS), RBDISTANCE, PVEC(3), RTEMP1(3,3), RTEMP2(3,3), SITEDIST
      DOUBLE PRECISION   :: QCUMUL(4), QBEST(4), QA(4), QB(4), QBINV(4), QTMP(4), QI(4)
      DOUBLE PRECISION   :: COORDSAS(3*NATOMS), COORDSBS(3*NATOMS), T(3*NATOMS)

      COORDSAS(1:3*NATOMS) = COORDSA(1:3*NATOMS) ! TO TRACE ERROR, SEE AT THE END
      COORDSBS(1:3*NATOMS) = COORDSB(1:3*NATOMS) ! TO TRACE ERROR, SEE AT THE END

      CMBX = 0.0D0; CMBY = 0.0D0; CMBZ = 0.0D0
      DO J1 = 1, NATOMS
         J2 = 3*J1
         CMBX = CMBX + COORDSB(J2-2)
         CMBY = CMBY + COORDSB(J2-1)
         CMBZ = CMBZ + COORDSB(J2)
       ENDDO
      CMBX = CMBX/NATOMS; CMBY = CMBY/NATOMS; CMBZ = CMBZ/NATOMS
!
!     BRING COORDSB INTO STANDARD ORIENTATION WITH RESPECT TO THE SITE POSITION
!     THE STANDARD ORIENTATION NEEDS TO BE DONE FOR THE SITES IF WE ARE GOING TO IDENTIFY
!     PERMUTATION-INVERSION ISOMERS WITH RESPECT TO THE SITES METRIC!
!
!     ----------------------------------------------------------------------------------------------
!     !
!     CALL POTENTIAL(COORDSB,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
!     PRINT '(2(A,F25.15))','BEFORE ENERGYB =',ENERGY,' RMS=',RMS
!     CALL POTENTIAL(DUMMYB,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
!     PRINT '(2(A,F25.15))',' AFTER ENERGYB =',ENERGY,' RMS=',RMS
!     !
!     ---------------------------------------------------------------------------------------------- 
!
      INVERT = 1

      DUMMYB(1:3*NATOMS) = COORDSB(1:3*NATOMS)
      DUMMYC(1:3*NATOMS) = DUMMYB(1:3*NATOMS)
      CALL ORIENTA(DUMMYC,DUMMY,NORBIT1,1,NORBIT2,1,NATOMS,QB,DEBUG)
      CALL QROTMAT(QB,ROTB)
      DUMMYB(1:3*NATOMS) = DUMMY(1:3*NATOMS)

      DBEST    = 1.0D100
60    NCHOOSE1 = 0
65    NCHOOSE1 = NCHOOSE1+1
40    NCHOOSE2 = 0
30    NCHOOSE2 = NCHOOSE2+1

      DUMMYA(1:3*NATOMS) = COORDSA(1:3*NATOMS)

      DO J1 = 1, NATOMS
         ALLPERM(J1) = J1
      ENDDO
 
!     THE OPTIMAL ALIGNMENT RETURNED BY MINPERDIST IS A LOCAL MINIMUM, BUT MAY NOT
!     BE THE GLOBAL MINIMUM. CALLING RBORIENT FIRST SHOULD PUT PERMUTATIONAL ISOMERS
!     INTO A STANDARD ALIGNMENT AND SPOT THE GLOBAL MINIMUM ZEDRO DISTANCE IN ONE
!     GO. HOWEVER, WE ALSO NEED TO CYCLE OVER EQUIVALENT ATOMS IN ORBITS USING NCHOOSE2.
!
!     PROBLEMS CAN OCCUR IF WE DON'T USE ALL THE ATOMS SPECIFIED BY NORBIT1 AND NORBIT2
!     BECAUSE OF THE NUMERICAL CUTOFFS EMPLOYED IN MYORIENT. WE COULD MISS THE
!     RIGHT ORIENTATION! 
!
!     IF WE USE RBORIENT TO PRODUCE PARTICULAR ORIENTATIONS THEN WE END UP ALIGNING 
!     COORDSA NOT WITH COORDSB BUT WITH THE STANDARD ORIENTATION OF COORDSB IN DUMMYB.
!     WE NOW DEAL WITH THIS BY TRACKING THE COMPLETE TRANSFORMATION, INCLUDING THE
!     CONTRIBUTION OF MYORIENT USING ROTB AND ROTINVB.
!
      DUMMYC(1:3*NATOMS) = INVERT*DUMMYA(1:3*NATOMS)

      CALL ORIENTA(DUMMYC,DUMMY,NORBIT1,NCHOOSE1,NORBIT2,NCHOOSE2,NATOMS,QA,DEBUG)
      CALL QROTMAT(QA,ROTA)
      DUMMYA(1:3*NATOMS)=DUMMY(1:3*NATOMS)

!      PRINT *, 'DUMMYA'
!      PRINT *, DUMMYA
!      PRINT *, 'DUMMYB'
!      PRINT *, DUMMYB
!      STOP
!     ----------------------------------------------------------------------------------------------
!     !
!      CALL POTENTIAL(DUMMYB,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
!      PRINT '(2(A,F25.15))', ' ENERGYBD =',ENERGY,' RMS=',RMS
!      CALL POTENTIAL(DUMMYA,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
!      PRINT '(2(A,F25.15))', ' ENERGYAD =',ENERGY,' RMS=',RMS
!
!      DISTANCE = 0.D0
!      DO J1 = 1, 3*NRB
!         DISTANCE = DISTANCE + (DUMMYA(J1) - DUMMYB(J1))**2
!      ENDDO
!     !
!     ----------------------------------------------------------------------------------------------
      
      DISTANCE = 0.D0
      DO J1 = 1, 3*NATOMS
         DISTANCE = DISTANCE + (DUMMYA(J1) - DUMMYB(J1))**2
      ENDDO

!     ----------------------------------------------------------------------------------------------
!     !
!      WRITE(975,'(I6)') NRB*NRBSITES
!      WRITE(975,'(A,2I6,A,F20.10)') 'A SITES BEFORE ROTATION'
!      WRITE(975,'(A,3F20.10)') ('LA ',XA(3*(J3-1)+1:3*(J3-1)+3),J3=1,NRB*NRBSITES)
!      WRITE(975,'(I6)') NRB*NRBSITES
!      WRITE(975,'(A,2I6,A,F20.10)') 'B SITES BEFORE ROTATION'
!      WRITE(975,'(A,3F20.10)') ('LA ',XB(3*(J3-1)+1:3*(J3-1)+3),J3=1,NRB*NRBSITES)
!     !
!     ----------------------------------------------------------------------------------------------

      IF (DEBUG) PRINT '(A,G20.10)',' MINPERMDISTRBCOM> AFTER INITIAL CALL TO RBSITESORIENT DISTANCE=',SQRT(DISTANCE)
!
!     BIPARTITE MATCHING ROUTINE FOR PERMUTATIONS. COORDINATES IN DUMMYB DO NOT CHANGE
!     BUT THE COORDINATES IN DUMMYA DO. DISTANCE IS THE DISTANCE IN THIS CASE.
!     WE RETURN TO LABEL 10 AFTER EVERY ROUND OF PERMUTATIONAL/ORIENTATIONAL ALIGNMENT
!     UNLESS WE HAVE CONVERGED TO THE IDENTITY PERMUTATION.
!
!     ATOMS ARE NOT ALLOWED TO APPEAR IN MORE THAN ONE GROUP.
!     THE MAXIMUM NUMBER OF PAIR EXCHANGES ASSOCIATED WITH A GROUP IS TWO.
!
      NTRIES = 0
!
!     QCUMUL IS A QUATERNION CONTAINING THE INFORMATION OF ACCUMULATED ROTATION THAT RELATES THE  
!     ORIGINAL DUMMYA OBTAINED FROM COORDSA TO THE FINAL ONE. RMATCUMUL IS THE CORRESPONDING
!     ROTATION MATRIX. INITIALIZE QCUMUL, SO AS RMATCUMUL.
!
      QCUMUL(1:4) = (/1.D0, 0.D0, 0.D0, 0.D0/)
      RMATCUMUL(:,:) = 0.D0; RMATCUMUL(1,1) = 1.D0; RMATCUMUL(2,2) = 1.D0; RMATCUMUL(3,3) = 1.D0
10    CONTINUE

      NTRIES = NTRIES + 1
      NDUMMY = 1

      DO J1 = 1, NATOMS
         NEWPERM(J1) = J1
      ENDDO

!     ALLPERM SAVES THE PERMUTATION FROM THE PREVIOUS CYCLE.
!     NEWPERM CONTAINS THE PERMUTATION FOR THIS CYCLE, RELATIVE TO THE IDENTITY.
!     SAVEPERM IS TEMPORARY STORAGE FOR NEWPERM.
!     NEWPERM MUST BE APPLIED TO ALLPERM AFTER THE LOOP OVER NPERMGROUP AND CORRESPONDING SWAPS.

      DO J1 = 1, NPERMGROUP

         PATOMS = NPERMSIZE(J1)
          IF (PATOMS.GT.NATOMS) THEN
             PRINT '(A,I6,A,I6,A,I6)',' MINPERMDISTRBCOM> ERROR *** NUMBER OF PERMUTABLE SITES IN GROUP ',J1,' IS ',PATOMS, &
  &                                   ' WHICH EXCEEDS THE NUMBER OF ATOMS ',NATOMS
             STOP
         ENDIF

         DO J2 = 1, PATOMS
            PDUMMYA(3*(J2-1)+1)=DUMMYA(3*(NEWPERM(PERMGROUP(NDUMMY+J2-1))-1)+1)
            PDUMMYA(3*(J2-1)+2)=DUMMYA(3*(NEWPERM(PERMGROUP(NDUMMY+J2-1))-1)+2)
            PDUMMYA(3*(J2-1)+3)=DUMMYA(3*(NEWPERM(PERMGROUP(NDUMMY+J2-1))-1)+3)
            PDUMMYB(3*(J2-1)+1)=DUMMYB(3*(NEWPERM(PERMGROUP(NDUMMY+J2-1))-1)+1)
            PDUMMYB(3*(J2-1)+2)=DUMMYB(3*(NEWPERM(PERMGROUP(NDUMMY+J2-1))-1)+2)
            PDUMMYB(3*(J2-1)+3)=DUMMYB(3*(NEWPERM(PERMGROUP(NDUMMY+J2-1))-1)+3)
         ENDDO
!
!     ALL PERMUTATIONS WITHIN THIS GROUP OF SIZE NPERMSIZE(J1) ARE NOW TRIED.
!     NOTE THAT WE ARE JUST USING A METRIC BASED ON THE RIGID BODY CENTRE OF MASS COORDINATES HERE!
!

         CALL MINPERM(PATOMS, PDUMMYB, PDUMMYA, BOXLX, BOXLY, BOXLZ, BULKT, LPERM, LDISTANCE, DIST2, WORSTRAD)

         SAVEPERM(1:NATOMS)=NEWPERM(1:NATOMS)
         DO J2=1,PATOMS
            SAVEPERM(PERMGROUP(NDUMMY+J2-1))=NEWPERM(PERMGROUP(NDUMMY+LPERM(J2)-1))
         ENDDO
!
! UPDATE PERMUTATION OF ASSOCIATED ATOMS, IF ANY.
! WE MUST DO THIS AS WE GO ALONG, BECAUSE THESE ATOMS COULD MOVE IN MORE THAN
! ONE PERMUTATIONAL GROUP NOW.
!
         IF (NSETS(J1) > 0) THEN
            DO J2 = 1, PATOMS
               DO J3 = 1, NSETS(J1)
                  SAVEPERM(SETS(PERMGROUP(NDUMMY+J2-1),J3))=SETS(NEWPERM(PERMGROUP(NDUMMY+LPERM(J2)-1)),J3)
               ENDDO
            ENDDO
         ENDIF
         NDUMMY = NDUMMY + NPERMSIZE(J1)
         NEWPERM(1:NATOMS) = SAVEPERM(1:NATOMS)

      ENDDO
!
!     DUE TO POSSIBLE SWAPS ABOVE, WE NEED TO UPDATE THE OVERALL PERMUTATION HERE.
!
!     UPDATE THE OVERALL PERMUTATION HERE.
!
      DO J1=1,NATOMS
         SAVEPERM(ALLPERM(J1))=ALLPERM(NEWPERM(J1))
      ENDDO

      ALLPERM(1:NATOMS) = NEWPERM(1:NATOMS)
      DUMMY(1:3*NATOMS) = DUMMYA(1:3*NATOMS)
      NPERM    = 0
      DISTANCE = 0.0D0
!
!     NOW PERMUTE THE ROTATIONAL COORDINATES TO CORRESPOND TO THE CENTRE OF MASS COORDINATE PERMUTATION.
!
!      DO J1 = (NATOMS/2)+1, NATOMS
!         IF (PERM(J1).NE.J1) THEN
!            PRINT '(A,2I8)',' MINPERMDISTRBCOM> ERROR - J1,PERM(J1) SHOULD BE EQUAL', J1,PERM(J1)
!         ENDIF
!         PERM(J1) = PERM(J1-(NATOMS/2)) + (NATOMS/2)
!      ENDDO

      DO J3 = 1, NATOMS

         DUMMYA(3*(J3-1)+1) = DUMMY(3*(NEWPERM(J3)-1)+1)
         DUMMYA(3*(J3-1)+2) = DUMMY(3*(NEWPERM(J3)-1)+2)
         DUMMYA(3*(J3-1)+3) = DUMMY(3*(NEWPERM(J3)-1)+3)

         IF (J3.NE.NEWPERM(J3)) THEN

            IF (DEBUG) WRITE(*,'(A,I5,A,I5)') ' MINPERMDISTRBCOM> MOVE POSITION ',NEWPERM(J3),' TO ',J3
            NPERM = NPERM + 1

         ENDIF

      ENDDO

      DO J1 = 1, 3*NATOMS
         DISTANCE = DISTANCE + (DUMMYA(J1) - DUMMYB(J1))**2
      ENDDO
!
!     FURTHER ALIGNMENT. COORDINATES IN DUMMYA ARE RESET BY RBMINDIST (SECOND ARGUMENT).
!     MUST ALLOW AT LEAST ONE CALL TO RBMINDIST IN CASE THE STANDARD ORIENTATION RESULT IS TERRIBLE
!     BUT GIVES ZERO PERMUTATIONS!
!     WE TRY INTERNAL SYMMETRY OPERATIONS FIRST FOR EACH RIGID BODY, THEN MINIMISE THE
!     DISTANCE FURTHER (IF POSSIBLE) .
!  
      IF ((NPERM.NE.0) .OR. (NTRIES.EQ.1)) THEN 
!
!     NOW IF RBSYMT IS .TRUE. WE SHOULD MINIMISE THE DISTANCE METRIC FOR ALL THE RIGID BODIES
!     BY CONSIDERING ALL THE ALLOWED SYMMETRY OPERATIONS FOR EACH RIGID BODY IN TURN.
!     DO THIS BY ALTERING THE ORIENTATIONAL COORDINATES IN DUMMYA FOR EACH RIGID BODY IN TURN,
!     AFTER SAVING THEM, CALCULATE THE NEW DISTANCE USING THE ALL SITES METRIC, AND ACCEPT THE
!     NEW ORIENTATION IF THE DISTANCE IS LOWER.
!     COULD GENERALISE TO MORE THAN ONE SORT OF RIGID BODY AS WELL.
! 
!     THIS CALL ALIGNS THE OVERALL ORIENTATION WITH RESPECT TO THE SITES METRIC.
!     INTERNAL SYMMETRY OPERATIONS OF THE RIGID BODIES ARE NOT CONSIDERED.
!
            CALL MINDISTA(DUMMYB,DUMMYA,NATOMS,DISTANCE,Q2,DEBUG)

            CALL QROTMAT(Q2,RMAT)

            DISTANCE  = DISTANCE*DISTANCE

!     ACCUMULATE ROTATION

         CALL QROTQ(Q2,QCUMUL)

         RMATCUMUL(:,:) = MATMUL(RMAT,RMATCUMUL)  ! ALSO IN THE FORM OF THE ROTATIONA MATRIX

         IF (NTRIES .LT. MAXIMUMTRIES) THEN
            GOTO 10
         ELSE ! PREVENT INFINITE LOOP
            IF (INVERT == -1) THEN
               IF(DEBUG) PRINT '(A)',' MINPERMDISTRBCOM> WARNING - NUMBER OF TRIES EXCEEDED, GIVING UP'
            ELSE
               PRINT '(A)',' MINPERMDISTRBCOM> WARNING - NUMBER OF TRIES EXCEEDED, GIVING UP'
            ENDIF
         ENDIF

      ENDIF

      IF (DISTANCE .LT. DBEST) THEN

         DBEST                = DISTANCE
         XBEST(1:3*NATOMS)    = DUMMYA(1:3*NATOMS)
         QTMP(1:4)  = QA(1:4)
         CALL QROTQ(QCUMUL,QTMP)
         QBEST(1:4) = QTMP(1:4)
         QBINV(1:4) = (/QB(1), -QB(2:4)/)
         RMATBEST(:,:) = MATMUL(RMATCUMUL,ROTA)
         IF (INVERT == -1) THEN
            RMATBEST(:,:) =-RMATBEST(:,:)
         ENDIF

      ENDIF
!
! IF GEOMDIFFTOLL IS SET TOO SLOPPY WE COULD MISS THE BEST SOLUTION BY EXITING VIA THE LINE
! BELOW WITHOUT TRYING OTHER ORBITS. TURN OFF THIS ESCAPE?!
! THE TURN OFF SEEMS TO BE A BUG FOR NON-RBPERM RUNS. LJ38 TESTS FAIL WITH ALL COMPILERS, PRODUCING
! A "DISTANCE IS ZERO: THIS SHOULD NOT HAPPEN" ERROR!!! DJW
! PUT THE ESCAPE BACK IN FOR NOW.
!
      IF (SQRT(DBEST).LT.GEOMDIFFTOL) GOTO 50
      IF (NCHOOSE2.LT.NORBIT2) GOTO 30
      IF (NCHOOSE1.LT.NORBIT1) GOTO 65
      IF (EFIELDT) GOTO 50
!      GOTO 50  !!!! THIS IS A WORKAROUND TO PREVENT PROBLEMS WITH MAKING THE INVERTED STRUCTURE. DJW

      IF ((NCHOOSE2.EQ.NORBIT2).AND.(NCHOOSE1.EQ.NORBIT1).AND.(INVERT.EQ.1)) THEN
!
!     DON'T TRY INVERSION FOR BULK OR CHARMM OR AMBER OR FROZEN ATOMS
!
         IF (DEBUG) PRINT '(A)',' MINPERMDISTRBCOM> INVERTING GEOMETRY FOR COMPARISON WITH TARGET'
         INVERT=-1
         GOTO 60

      ENDIF

!
50    DISTANCE = DBEST       ! SQUARED DISTANCE
!
!     XBEST CONTAINS THE BEST ALIGNMENT OF A COORDINATES FOR THE ORIENTATION OF B COORDINATES IN DUMMYB.
!     ROTATE XBEST BY ROTINVBBEST TO PUT IN BEST CORRESPONDENCE WITH COORDSB, UNDOING THE REORIENTATION TO  
!     DUMMYB FROM RBORIENT. 
!     WE SHOULD GET THE SAME RESULT FOR ROTINVBBEST * RMATBEST * (COORDSA-CMA), 
!     WHERE RMATBEST = +/- RMATCUMUL * ROTA FOR THE BEST ALIGNMENT 
!     (ASIDE FROM A POSSIBLE PERMUTATION OF THE ATOM ORDERING)
!

!      CALL QROTMAT(QBINV, ROTBINV)

      ROTBINV = TRANSPOSE(ROTB)
 
      DO J1 = 1, NATOMS
         J2 = 3*J1
         XBEST(J2-2:J2) = MATMUL(ROTBINV,XBEST(J2-2:J2))
      ENDDO

      DO J1 = 1, NATOMS
         XBEST(3*(J1-1)+1) = XBEST(3*(J1-1)+1) + CMBX
         XBEST(3*(J1-1)+2) = XBEST(3*(J1-1)+2) + CMBY
         XBEST(3*(J1-1)+3) = XBEST(3*(J1-1)+3) + CMBZ
      ENDDO

      XDUMMY = 0.D0
      DO J1 = 1, 3*NATOMS
         XDUMMY = XDUMMY + (XBEST(J1) - COORDSB(J1))**2
      ENDDO

      IF (ABS(SQRT(XDUMMY)-SQRT(DISTANCE)).GT.GEOMDIFFTOL) THEN
          PRINT '(2(A,F20.10))','MINPERMDISTRBCOM> ERROR *** DISTANCE BETWEEN TRANSFORMED XBEST AND COORDSB=',  &
     &    SQRT(XDUMMY),  ' SHOULD BE ', SQRT(DISTANCE)
          PRINT *, 'COORDSBCOM'
          DO J1 = 1, NATOMS
             J2 = 3*J1
             PRINT *, COORDSBS(J2-2), COORDSBS(J2-1), COORDSBS(J2)
          ENDDO 
          PRINT *, 'COORDSACOM'
          DO J1 = 1, NATOMS
             J2 = 3*J1
             PRINT *, COORDSAS(J2-2), COORDSAS(J2-1), COORDSAS(J2)
          ENDDO 

          STOP

      ENDIF

      CALL QROTQ(QBINV,QBEST)

      RMATBEST(:,:) = MATMUL(ROTBINV,RMATBEST)

      COORDSA(1:3*NATOMS)=XBEST(1:3*NATOMS) ! FINALLY, BEST COORDSA SHOULD INCLUDE PERMUTATIONS FOR DNEB INPUT!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!DEBUG
!      CALL POTENTIAL(COORDSA,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
!      PRINT '(2(A,F25.15))',' FINALLY A ENERGY=',ENERGY,' RMS=',RMS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!DEBUG

      END SUBROUTINE MINPERMDISTRBCOM

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE ORIENTA(X, T1S, NORBIT1, NCHOOSE1, NORBIT2, NCHOOSE2, NATOMS, Q2, DEBUG)

!     THIS SUBROUTINE PUTS THE CONFIGURATION, X, OF AN ATOMIC  SYSTEM INTO A STANDARD ALIGNMENT, T1.
!
      USE KEY, ONLY: EFIELDT

      IMPLICIT NONE
      INTEGER          :: NATOMS, I, J, J1, J2, JMAX1, JMAX2, NORBIT1, NCHOOSE1, NORBIT2, NCHOOSE2
      DOUBLE PRECISION :: X(3*NATOMS), T1(3*NATOMS)
      DOUBLE PRECISION :: XS(3*NATOMS), T1S(3*NATOMS), T2S(3*NATOMS), DIST(NATOMS)
      DOUBLE PRECISION :: AX(3), P(3), Q2(4), ROTM(3,3), ROTMINV(3,3)
      DOUBLE PRECISION :: THETA, THETAH, COST, SINT, COSTH, SINTH, ST, FCT
      DOUBLE PRECISION :: CMX, CMY, CMZ, DMAX, DUMMY, PROJ, DMAX2, CUTOFF1, DTEMP
      LOGICAL          :: DEBUG

      CUTOFF1 = 1.D-03
!
!     MOVE CENTRE OF MASS TO THE ORIGIN.
!
      CMX = 0.D0; CMY = 0.D0; CMZ = 0.D0
      DO I = 1, NATOMS
         J = 3*I
         CMX = CMX + X(J-2)
         CMY = CMY + X(J-1)
         CMZ = CMZ + X(J)
      ENDDO

      CMX = CMX/NATOMS; CMY = CMY/NATOMS; CMZ = CMZ/NATOMS
      DO I = 1, NATOMS
         J = 3*I
         XS(J-2) = X(J-2) - CMX
         XS(J-1) = X(J-1) - CMY
         XS(J)   = X(J)   - CMZ
      ENDDO

      DMAX    = -1.D0
      NORBIT1 = 1

      IF (EFIELDT) THEN
         T1S(:) = XS(:)
         GOTO 100
      ENDIF

!
!     FIND THE ATOM WHICH IS AT THE LARGEST DISTANCE FROM THE CENTRE OF MASS
!
      DO J1 = 1, NATOMS

         J = 3*J1
         DIST(J1) = SQRT(XS(J-2)**2 + XS(J-1)**2 + XS(J)**2)

         IF (ABS(DIST(J1) - DMAX) < CUTOFF1) THEN

            NORBIT1 = NORBIT1 + 1
            IF (NORBIT1 == NCHOOSE1) JMAX1 = J1

         ELSE IF (DIST(J1) > DMAX) THEN

            DMAX = DIST(J1)
            NORBIT1 = 1
            JMAX1 = J1

         ENDIF

      ENDDO

!     FOR TAGGED ATOMS, THE CHOICE OF THE FIRST ATOM MATTERS IF IT BELONGS TO AN ORBIT OF SIZE > 1.
!
      IF ((ABS(XS(3*JMAX1-2)) < 1.D-08) .AND. (ABS(XS(3*JMAX1-1)) < 1.D-08)) THEN

!
!     ATOM JMAX1 IS ALREADY ON THE Z AXIS!
!
         IF (XS(3*(JMAX1-1)+3) > 0.D0) THEN

            T1S(1:3*NATOMS) = XS(1:3*NATOMS)
            Q2(1:4)   = (/1.D0, 0.D0, 0.D0, 0.D0/)   ! IDENTITY OPERATION       
                  
         ELSE  ! ROTATE ABOUT THE X-AXIS BY \PI, DO NOT INVERT!

            Q2(1:4) = (/0.D0, 1.D0, 0.D0, 0.D0/)      ! THE CORRESPONDING QUATERNION

            DO J1 = 1, NATOMS
               J       = 3*J1
               T1S(J-2) = XS(J-2)
               T1S(J-1) =-XS(J-1)
               T1S(J)   =-XS(J)
            ENDDO

         ENDIF

      ELSE

!
!     ROTATE ALL ATOMS THROUGH AN ANGLE THETA ABOUT THE AXIS P(:), SO AS TO ROTATE THE RBSITE JMAX1 ONTO THE Z AXIS 
!
!       THETA = DACOS(XS(3*JMAX1)/DMAX)
!     FOR SLOPPY CUTOFFS WE CANNOT ASSUME THAT DMAX IS EXACTLY THE SAME FOR MEMBERS OF THE SAME ORBIT!

        THETA = DACOS(XS(3*JMAX1)/DIST(JMAX1))

!     THE AXIS IS ON THE XY-PLANE AND PERPENDICULAR TO THE PROJECTION OF THE TRANSLATION VECTOR OF THE RB JMAX1
!     ON THE XY-PLANE. THE AXIS IS SO CHOSEN THAT A RIGHT-HANDED ROTATION IN THE RIGHT-HANDED COORDINATE SYSTEM
!     WILL RESULT IN THE REQUIRED TRANSFORMATION.

         AX(3) = 0.D0
         FCT   = DSQRT(XS(3*JMAX1-2)**2 + XS(3*JMAX1-1)**2)
         AX(1) = THETA*XS(3*JMAX1-1)/FCT
         AX(2) =-THETA*XS(3*JMAX1-2)/FCT

         CALL ROTMAT(AX(:), ROTM(:,:))

!         CALL INVMTRX(ROTM, ROTMINV)

         THETAH  = 0.5D0*THETA
         AX(1:3) = AX(1:3)/SQRT(DOT_PRODUCT(AX(1:3),AX(1:3)))
         Q2(1:4) = (/COS(THETAH), SIN(THETAH)*AX(1:3)/)          ! THE CORRESPONDING QUATERNION 
          
         DO J1 = 1, NATOMS
            J2 = 3*J1
            T1S(J2-2:J2) = MATMUL(ROTM(:,:),XS(J2-2:J2))
         ENDDO   

      ENDIF
!
!     NOW FIND THE ATOM WITH THE LARGEST DISTANCE FROM THE Z-AXIS 
!
100   DMAX = -1.0D0

      DO J1 = 1, NATOMS
         J2 = 3*J1 
         DIST(J1) = SQRT(T1S(J2-2)**2 + T1S(J2-1)**2)
         IF (DIST(J1) > DMAX) DMAX = DIST(J1)
      ENDDO

      DMAX2 = -1.0D100
!
!     PROJ IS THE SUM OF THE X COMPONENTS. USE T2S AS A DUMMY IN ORDER NOT TO 
!     CHANGE T1 UNTIL WE HAVE DECIDED WHICH ATOM TO PUT IN THE XZ PLANE.
!
      DO J1 = 1, NATOMS

         IF (ABS(DIST(J1) - DMAX) < CUTOFF1) THEN

            T2S(1:3*NATOMS) = T1S(1:3*NATOMS)

            CALL ROTATMXZ(NATOMS, J1, T2S, PROJ, DIST, Q2, .FALSE.)

            IF (ABS(PROJ - DMAX2) < CUTOFF1) THEN
               NORBIT2 = NORBIT2+1
               IF (NORBIT2 == NCHOOSE2) THEN
                  JMAX2 = J1
                  DTEMP = PROJ
               ENDIF
            ELSE IF (PROJ > DMAX2) THEN
               NORBIT2 = 1
               DMAX2   = PROJ
               DTEMP   = PROJ
               JMAX2   = J1
            ENDIF

         ENDIF

      ENDDO

!
!     AND NOW ROTATE IT INTO THE XZ PLANE.
!

      CALL ROTATMXZ(NATOMS, JMAX2, T1S, DTEMP, DIST, Q2, .TRUE.)

      END SUBROUTINE ORIENTA

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE ROTATMXZ(NATOMS, JDO, T1S, PROJ, DIST, Q1, ROTT)

      USE KEY, ONLY : EFIELDT

      IMPLICIT NONE
      INTEGER          :: NATOMS, JDO, I, J, J2
      DOUBLE PRECISION :: T1S(3*NATOMS), PROJ, DIST(NATOMS), THETA, THETAH, COST, SINT, ST
      DOUBLE PRECISION :: COSTH, SINTH, FCT, P(3), Q2(4), Q1(4), RM(3,3), ROTM(3,3) 
      LOGICAL          :: ROTT

      J2     = 3*JDO

      IF (ABS(T1S(J2-1)) < 1.0D-8) THEN            ! ALREADY ON THE XZ PLANE

         Q2(1:4) = (/1.D0, 0.D0, 0.D0, 0.D0/) 

         IF (T1S(J2-2) < 0.D0) THEN                ! ROTATE ABOUT THE Z AXIS BY \PI, DO NOT INVERT!!

            Q2(1:4) = (/0.D0, 0.D0, 0.D0, 1.D0/) 

            DO I = 1, NATOMS
               J = 3*I
               T1S(J-2) =-T1S(J-2)
               T1S(J-1) =-T1S(J-1)
            ENDDO

         ENDIF

      ELSE

         COST   = T1S(J2-2)/DIST(JDO)       
         SINT   = T1S(J2-1)/DIST(JDO)
         THETA  =-ATAN2(SINT,COST)              ! THE NEGATIVE SIGN APPEARS AS THE ROTATION IS REVERSED  
         THETAH = 0.5D0*THETA
         Q2(1:4) = (/COS(THETAH), 0.D0, 0.D0, SIN(THETAH)/)     ! THE QUATERNION CORRESPONDING TO R_{Z}(-\THETA)
         RM(1:3,1:3) = 0.D0                           ! THE ROTATION MATRIX FOR R_{Z}(-\THETA)
         RM(1,1) = COST
         RM(2,2) = COST
         RM(3,3) = 1.D0
         RM(1,2) = SINT
         RM(2,1) =-SINT

         CALL QROTMAT(Q2,RM)

         DO I = 1, NATOMS
            IF (DIST(I) /= 0.0D0) THEN
               J = 3*I
               T1S(J-2:J) = MATMUL(RM(:,:),T1S(J-2:J))
            ENDIF  
         ENDDO

      ENDIF

      IF (ROTT) THEN

         IF (EFIELDT) THEN
            Q1(1:4) = Q2(1:4)
         ELSE
            CALL QROTQ(Q2,Q1)
         ENDIF

      ELSE

         PROJ = 0.D0

         DO I = 1, NATOMS 

            J = 3*I
            IF (T1S(J) > 1.0D-2) PROJ = PROJ + T1S(J-2)

         ENDDO

      ENDIF

      END SUBROUTINE ROTATMXZ

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE ROTATM(T, X, Q2, NATOMS)
!     TAKES THE SET OF COORDINATES T FOR NATOMS NUMBER OF ATOMS AND RETURNS X AFTER ROTATION VIA THE 
!     QUATERNION Q2 ABOUT THE ORIGIN 

      IMPLICIT NONE

      INTEGER          :: I, J, NATOMS
      DOUBLE PRECISION :: T(3*NATOMS), X(3*NATOMS), T1(1:3), Q2(4), RM(3,3) 
      DOUBLE PRECISION :: CMX, CMY, CMZ
!
!     MOVE CENTRE OF MASS TO THE ORIGIN.
!
      CMX = 0.D0; CMY = 0.D0; CMZ = 0.D0
      DO I = 1, NATOMS
         J = 3*I
         CMX = CMX + T(J-2)
         CMY = CMY + T(J-1)
         CMZ = CMZ + T(J)
      ENDDO
      CMX = CMX/NATOMS; CMY = CMY/NATOMS; CMZ = CMZ/NATOMS
      DO I = 1, NATOMS
         J      = 3*I
         X(J-2) = T(J-2) - CMX
         X(J-1) = T(J-1) - CMY
         X(J)   = T(J)   - CMZ
      ENDDO

!     EXTRACT THE ROTATION MATRIX CORRESPONDING TO ROTATION VIA Q

      CALL QROTMAT(Q2,RM)

      DO I = 1, NATOMS

         J        = 3*I
         T1(1:3)  = MATMUL(RM(:,:), X(J-2:J))
         X(J-2:J) = T1(1:3)  

      ENDDO

      END SUBROUTINE ROTATM

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE MINDISTA(RA,RB,NATOMS,DIST,RM,DEBUG)

!     RETURNS SQUARED DISTANCE DIST

      USE KEY, ONLY: EFIELDT

      IMPLICIT NONE

      INTEGER          :: J1, J2, J3, J4, NATOMS, NSIZE, JMIN, INFO
      DOUBLE PRECISION :: RA(3*NATOMS), RB(3*NATOMS), DIST, QMAT(4,4), TEMPA(9*NATOMS), XM, YM, ZM, XP, YP, ZP
      DOUBLE PRECISION :: DIAG(4), MINV, Q2(4), CMXA, CMYA, CMZA, CMXB, CMYB, CMZB
      DOUBLE PRECISION :: R(3), P(3), RM(3,3)
      DOUBLE PRECISION, ALLOCATABLE :: XA(:), XB(:)
      DOUBLE PRECISION :: ENERGY, VNEW(3*NATOMS), RMS, DUMMY
      LOGICAL          :: BULKT, PRESERVET, DEBUG

      NSIZE = NATOMS
      ALLOCATE(XA(3*NSIZE),XB(3*NSIZE))
      XA(1:3*NSIZE) = RA(1:3*NSIZE); XB(1:3*NSIZE) = RB(1:3*NSIZE)

!     MOVE CENTRE OF COORDINATES OF XA AND XB TO THE ORIGIN

      CMXA = 0.0D0; CMYA = 0.0D0; CMZA = 0.0D0
      DO J1 = 1, NSIZE
         J2 = 3*(J1-1)
         CMXA = CMXA + XA(J2+1)
         CMYA = CMYA + XA(J2+2)
         CMZA = CMZA + XA(J2+3)
      ENDDO
      CMXA = CMXA/NSIZE; CMYA = CMYA/NSIZE; CMZA = CMZA/NSIZE
      DO J1 = 1, NSIZE
         J2 = 3*(J1-1)
         XA(J2+1) = XA(J2+1) - CMXA
         XA(J2+2) = XA(J2+2) - CMYA
         XA(J2+3) = XA(J2+3) - CMZA
      ENDDO

      CMXB = 0.0D0; CMYB = 0.0D0; CMZB = 0.0D0
      DO J1 = 1, NSIZE
         J2 = 3*(J1-1)
         CMXB = CMXB + XB(J2+1)
         CMYB = CMYB + XB(J2+2)
         CMZB = CMZB + XB(J2+3)
      ENDDO
      CMXB = CMXB/NSIZE; CMYB = CMYB/NSIZE; CMZB = CMZB/NSIZE
      DO J1 = 1, NSIZE
         J2 = 3*(J1-1)
         XB(J2+1) = XB(J2+1) - CMXB
         XB(J2+2) = XB(J2+2) - CMYB
         XB(J2+3) = XB(J2+3) - CMZB
      ENDDO

      QMAT(1:4,1:4) = 0.0D0

      DO J1 = 1, NSIZE
         J2 = 3*(J1-1)
         XM = XA(J2+1) - XB(J2+1)
         YM = XA(J2+2) - XB(J2+2)
         ZM = XA(J2+3) - XB(J2+3)
         XP = XA(J2+1) + XB(J2+1)
         YP = XA(J2+2) + XB(J2+2)
         ZP = XA(J2+3) + XB(J2+3)

         IF (EFIELDT) THEN
            QMAT(1,1) = QMAT(1,1) + XM**2 + YM**2 + ZM**2
            QMAT(1,2) = QMAT(1,2) - XP*YM + XM*YP
            QMAT(2,2) = QMAT(2,2) + XP**2 + YP**2 + ZM**2
         ELSE
            QMAT(1,1) = QMAT(1,1) + XM**2 + YM**2 + ZM**2
            QMAT(1,2) = QMAT(1,2) - YP*ZM + YM*ZP
            QMAT(1,3) = QMAT(1,3) - XM*ZP + XP*ZM
            QMAT(1,4) = QMAT(1,4) - XP*YM + XM*YP
            QMAT(2,2) = QMAT(2,2) + YP**2 + ZP**2 + XM**2
            QMAT(2,3) = QMAT(2,3) + XM*YM - XP*YP
            QMAT(2,4) = QMAT(2,4) + XM*ZM - XP*ZP
            QMAT(3,3) = QMAT(3,3) + XP**2 + ZP**2 + YM**2
            QMAT(3,4) = QMAT(3,4) + YM*ZM - YP*ZP
            QMAT(4,4) = QMAT(4,4) + XP**2 + YP**2 + ZM**2
         ENDIF
      ENDDO

      IF (EFIELDT) THEN

!     QMAT IS SYMMETRIC; QMAT(2,1) = QMAT(1,2)

         MINV = 0.5D0*(QMAT(1,1) + QMAT(2,2) - SQRT(4.D0*QMAT(1,2)*QMAT(1,2) + (QMAT(1,1) - QMAT(2,2))**2.D0))
         Q2(1) = SQRT((MINV-QMAT(2,2))**2.D0/(QMAT(1,2)*QMAT(1,2) + (MINV-QMAT(2,2))**2.D0))
         Q2(2) = 0.D0
         Q2(3) = 0.D0
         Q2(4) = QMAT(1,2)*Q2(1)/(MINV - QMAT(2,2))

         IF (MINV < 0.0D0) THEN
            IF (ABS(MINV)< 1.0D-6) THEN
               MINV = 0.0D0
            ELSE
               PRINT '(A,G20.10,A)','NEWMINDIST> WARNING MINV IS ',MINV,' CHANGE TO ABSOLUTE VALUE'
               MINV = -MINV
            ENDIF
         ENDIF
       
      ELSE

        QMAT(2,1) = QMAT(1,2); QMAT(3,1) = QMAT(1,3); QMAT(3,2) = QMAT(2,3); QMAT(4,1) = QMAT(1,4)
        QMAT(4,2) = QMAT(2,4); QMAT(4,3) = QMAT(3,4)
        CALL DSYEV('V','U',4,QMAT,4,DIAG,TEMPA,9*NATOMS,INFO)

        IF (INFO /= 0) PRINT '(A,I6,A)','MINDISTA> WARNING - INFO=',INFO,' IN DSYEV'

        MINV = 1.0D100
        DO J1 = 1,4
           IF (DIAG(J1).LT.MINV) THEN
              JMIN = J1
              MINV = DIAG(J1)
           ENDIF
        ENDDO
        IF (MINV < 0.0D0) THEN
           IF (ABS(MINV)< 1.0D-6) THEN
              MINV = 0.0D0
           ELSE
              PRINT '(A,G20.10,A)','MINDISTA> WARNING MINV IS ',MINV,' CHANGE TO ABSOLUTE VALUE'
              MINV = -MINV
           ENDIF
        ENDIF

        Q2(1) = QMAT(1,JMIN); Q2(2) = QMAT(2,JMIN); Q2(3) = QMAT(3,JMIN); Q2(4) = QMAT(4,JMIN)

      ENDIF

      DIST = MINV

      DO J1 = 1, NATOMS
         J2 = 3*(J1-1)
         RB(J2+1) = RB(J2+1) - CMXB
         RB(J2+2) = RB(J2+2) - CMYB
         RB(J2+3) = RB(J2+3) - CMZB
      ENDDO

      CALL NEWROTGEOMA(NATOMS,RB,Q2,RM,CMXA,CMYA,CMZA)

      DEALLOCATE(XA,XB)

      END SUBROUTINE MINDISTA

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE NEWROTGEOMA(NATOMS,COORDS,Q2,RM,CX,CY,CZ)

      IMPLICIT NONE

      INTEGER          :: I, J, NATOMS
      DOUBLE PRECISION :: COORDS(3*NATOMS), RM(3,3), CX, CY, CZ, R(3), P(3), Q1(4), Q2(4), Q(4)
      DOUBLE PRECISION :: THETA, THETAH, ST, FCT

!     RMAT CONTAINS THE MATRIX THAT MAPS RB ONTO THE BEST CORRESPONDENCE WITH RA

      CALL QROTMAT(Q2,RM)

      DO I = 1, NATOMS

         J    = 3*(I-1)
         R(:) = MATMUL(RM(:,:), COORDS(J+1:J+3))

         COORDS(J+1) = R(1) + CX
         COORDS(J+2) = R(2) + CY
         COORDS(J+3) = R(3) + CZ

      ENDDO

      END SUBROUTINE NEWROTGEOMA
