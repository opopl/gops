      SUBROUTINE FINALIO
      USE commons
      use modamber
      use modamber9, only : coords1,lcrd,ih,m04,natom,ambfinalio_node
      use pymodule, only : SITECOORDS,ELLST1,ELLMAT
      use qmodule
      USE modcharmm
      USE AMHGLOBALS, ONLY:NMRES,IRES

      IMPLICIT NONE

      INTEGER III, I3,  GLY_COUNT, ID, NUMCRD, NUMPRO, NCPHST
      INTEGER J1, J2, J3, J4, J5, NTYPEA, MYUNIT2, I1, NDUMMY, MYUNIT3, NC, NRBS1, NRBS2
      DOUBLE PRECISION EPSAB, EPSBB, SIGAB, SIGBB, RBCOORDS(NRBSITES*3), DCOORDS(3*NATOMS)
      DOUBLE PRECISION P3(3,3), P(3), DU(3), RMI(3,3), DRMI(3,3), PI, PHI, THT, CHI
      DOUBLE PRECISION, ALLOCATABLE :: XCOORDS(:), YCOORDS(:)
      DOUBLE PRECISION EulerPhi,EulerPsi,EulerTheta,EulerThetadeg,EulerPhiDeg,EulerPsiDeg  ! Euler angles for ellipsoids of revolution

      DOUBLE PRECISION EPS2, RAD, HEIGHT,sumx,sumy,sumz
      LOGICAL :: GTEST


      INTEGER COUNTT
      DOUBLE PRECISION  PPPCORD(NMRES*3*3,3,3,5)
      EXTERNAL NUM_TO_CHAR
      PI = 4.D0*DATAN(1.D0)

      NUMPRO = 1
      NUMCRD = 3

      ALLOCATE(DBNAME(NSAVE))

     
      DO J1=1,NSAVE
         WRITE(DBNUM,*) J1
         DBNAME(J1)='dbase.'//TRIM(ADJUSTL(DBNUM))
      ENDDO
      DCOORDS(1:3*NATOMS)=0.0D0

       IF (AMHT) THEN
         OPEN(UNIT=26,FILE='movie_gmin',STATUS='UNKNOWN')
         WRITE(26,334)NMRES,NUMCRD,NUMPRO,NSAVE
334      FORMAT(4(I8,1X),' NMRES NMCRD NUMPRO NMSNAP')
       ENDIF


      IF (MPIT) THEN
         WRITE (ISTR, '(I10)') MYUNIT-22980+1 
         MYUNIT2=(MYUNIT-22980+1)+100
         MYFILENAME2="lowest."//trim(adjustl(istr))
         OPEN(MYUNIT2,FILE=trim(adjustl(MYFILENAME2)), STATUS="unknown", form="formatted")
         IF (CHRMMT) THEN
            DO J1=1,NSAVE
               WRITE(DBNUM,*) J1
               DBNAME(J1)='dbase.'//TRIM(ADJUSTL(ISTR))//'.'//TRIM(ADJUSTL(DBNUM))             
            ENDDO
         ENDIF
         IF (CSMT.AND.(.NOT.SYMMETRIZECSM)) THEN
            MYUNIT3=(MYUNIT+22980+2)+111
            MYFILENAME3="CSMav."//trim(adjustl(istr))
         ENDIF
      ELSE
         MYUNIT2=25 
         OPEN(MYUNIT2,FILE='lowest',STATUS='UNKNOWN')
         IF (CSMT.AND.(.NOT.SYMMETRIZECSM)) THEN
            MYUNIT3=26 
            OPEN(MYUNIT3,FILE='CSMav.xyz',STATUS='UNKNOWN')
         ENDIF
      ENDIF

      DO J1=1,NSAVE
         IF (AMHT) THEN
            OPEN(UNIT=27,FILE='movie_gmin.'//COUNTTT//'.pdb',STATUS='UNKNOWN')
         ENDIF

         IF (RGCL2.OR.ARNO) THEN 
            WRITE(MYUNIT2,*) NATOMS+2
         ELSE IF (ELLIPSOIDT.OR.LJCAPSIDT.OR.GBT.OR.GBDT) THEN
            WRITE(MYUNIT2,*) NATOMS/2
         ELSE IF (AMHT) THEN
            WRITE(MYUNIT2,*) NMRES*3
         ELSE
            WRITE(MYUNIT2,*) NATOMS
         ENDIF
         WRITE(MYUNIT2,10) J1, QMIN(J1), FF(J1)
10       FORMAT('Energy of minimum ',I6,'=',F20.10,' first found at step ',I8)
         IF (MSORIGT.OR.FRAUSIT) THEN
            WRITE(MYUNIT2,20) (QMINP(J1,J2),J2=1,3*(NATOMS-NS))
20          FORMAT('Si',3F20.10)
         ELSE IF (MSTRANST) THEN
            WRITE(MYUNIT2,20) (QMINP(J1,J2),J2=1,3*(NATOMS-NS))
         ELSE IF (RGCL2) THEN
            WRITE(MYUNIT2,'(A,F20.10)') 'Cl 0.0 0.0 ', 0.995D0
            WRITE(MYUNIT2,'(A,F20.10)') 'Cl 0.0 0.0 ',-0.995D0
            WRITE(MYUNIT2,60) (QMINP(J1,J2),J2=1,3*(NATOMS-NS))
60          FORMAT('AR ',3F20.10)
         ELSE IF (AMHT) THEN
            WRITE(26,683)NUMCRD,j1,NUMCRD,REAL(NUMCRD),NUMCRD
683         FORMAT(3(I6,1X),F8.4,1X,I5,' STUCT SNAP T T TID')

            GLY_COUNT = 0

            DO 1964 III = 1,NMRES
               IF (IRES(III).EQ.8) THEN
                  PPPCORD(III, 1, 1, 1) = REAL(QMINP(J1,9*(III-1)+1- GLY_COUNT*3)) !  CA X
                  PPPCORD(III, 2, 1, 1) = REAL(QMINP(J1,9*(III-1)+2- GLY_COUNT*3)) !  CA Y
                  PPPCORD(III, 3, 1, 1) = REAL(QMINP(J1,9*(III-1)+3- GLY_COUNT*3)) !  CA Z
                  PPPCORD(III, 1, 1, 2) = REAL(QMINP(j1,9*(III-1)+1- GLY_COUNT*3)) !  CB X
                  PPPCORD(III, 2, 1, 2) = REAL(QMINP(J1,9*(III-1)+2- GLY_COUNT*3)) !  CB Y
                  PPPCORD(III, 3, 1, 2) = REAL(QMINP(J1,9*(III-1)+3- GLY_COUNT*3)) !  CB Z
                  PPPCORD(III, 1, 1, 3) = REAL(QMINP(J1,9*(III-1)+4- GLY_COUNT*3)) !  O X
                  PPPCORD(III, 2, 1, 3) = REAL(QMINP(J1,9*(III-1)+5- GLY_COUNT*3)) !  O Y
                  PPPCORD(III, 3, 1, 3) = REAL(QMINP(J1,9*(III-1)+6- GLY_COUNT*3)) !  O Z

                  WRITE(MYUNIT2,31) QMINP(J1,9*(III-1)+1-GLY_COUNT*3),QMINP(J1,9*(III-1)+2-GLY_COUNT*3),
     &                                                                QMINP(J1,9*(III-1)+3-GLY_COUNT*3)
                  WRITE(MYUNIT2,31) QMINP(J1,9*(III-1)+1-GLY_COUNT*3),QMINP(J1,9*(III-1)+2-GLY_COUNT*3),
     &                                                                QMINP(J1,9*(III-1)+3-GLY_COUNT*3)
                  WRITE(MYUNIT2,31) QMINP(J1,9*(III-1)+4-GLY_COUNT*3),QMINP(J1,9*(III-1)+5-GLY_COUNT*3),
     &                                                                QMINP(J1,9*(III-1)+6-GLY_COUNT*3)
31                FORMAT('AM',3G25.15)

                  GLY_COUNT = GLY_COUNT +1
               ELSE
                  PPPCORD(III, 1, 1, 1) = REAL(QMINP(J1,9*(iii-1)+1- GLY_COUNT*3)) !  CA X
                  PPPCORD(III, 2, 1, 1) = REAL(QMINP(J1,9*(III-1)+2- GLY_COUNT*3)) !  CA Y
                  PPPCORD(III, 3, 1, 1) = REAL(QMINP(J1,9*(III-1)+3- GLY_COUNT*3)) !  CA Z
                  PPPCORD(III, 1, 1, 2) = REAL(QMINP(J1,9*(III-1)+4- GLY_COUNT*3)) !  CB X
                  PPPCORD(III, 2, 1, 2) = REAL(QMINP(J1,9*(III-1)+5- GLY_COUNT*3)) !  CB Y
                  PPPCORD(III, 3, 1, 2) = REAL(QMINP(J1,9*(III-1)+6- GLY_COUNT*3)) !  CB Z
                  PPPCORD(III, 1, 1, 3) = REAL(QMINP(J1,9*(III-1)+7- GLY_COUNT*3)) !  O X
                  PPPCORD(III, 2, 1, 3) = REAL(QMINP(J1,9*(III-1)+8- GLY_COUNT*3)) !  O Y
                  PPPCORD(III, 3, 1, 3) = REAL(QMINP(J1,9*(III-1)+9- GLY_COUNT*3)) !  O Z
 
                  WRITE(MYUNIT2,31) QMINP(J1,9*(III-1)+1-GLY_COUNT*3),QMINP(J1,9*(III-1)+2-GLY_COUNT*3),
     &                                                                QMINP(J1,9*(III-1)+3-GLY_COUNT*3)
                  WRITE(MYUNIT2,31) QMINP(J1,9*(III-1)+4-GLY_COUNT*3),QMINP(J1,9*(III-1)+5-GLY_COUNT*3),
     &                                                                QMINP(J1,9*(III-1)+6-GLY_COUNT*3)
                  WRITE(MYUNIT2,31) QMINP(J1,9*(III-1)+7-GLY_COUNT*3),QMINP(J1,9*(III-1)+8-GLY_COUNT*3),
     &                                                                QMINP(J1,9*(III-1)+9-GLY_COUNT*3)
               ENDIF
1964        CONTINUE

            DO 526 III=1,NMRES
               WRITE(26,632)(PPPCORD(III,I3,1,1),I3=1,3),(PPPCORD(III,I3,1,2),I3=1,3),(PPPCORD(III,I3,1,3),I3=1,3)
632            FORMAT('CA: ',3(F8.3,1X),'CB: ',3(F8.3,1X),'OX: ', 3(F8.3,1X))
526         CONTINUE

            DO III = 1+1, NMRES
               PPPCORD(III,1,1,4) =
     &          0.4831806D0*PPPCORD(III-1,1,1,1) + 0.7032820D0*PPPCORD(III,1,1,1) - 0.1864626D0*PPPCORD(III-1,1,1,3)
               PPPCORD(III,2,1,4) =
     &          0.4831806D0*PPPCORD(III-1,2,1,1) + 0.7032820D0*PPPCORD(III,2,1,1) - 0.1864626D0*PPPCORD(III-1,2,1,3)
               PPPCORD(III,3,1,4) =
     &          0.4831806D0*PPPCORD(III-1,3,1,1) + 0.7032820d0*PPPCORD(III,3,1,1) - 0.1864626D0*PPPCORD(III-1,3,1,3)
            ENDDO

            DO III = 1, NMRES-1
               PPPCORD(III,1,1,5) =
     &          0.4436538d0*PPPCORD(III,1,1,1)+0.2352006D0*PPPCORD(III+1,1,1,1)+0.3211455D0*PPPCORD(III,1,1,3)
               PPPCORD(III,2,1,5) =
     &          0.4436538D0*PPPCORD(III,2,1,1)+0.2352006D0*PPPCORD(III+1,2,1,1)+0.3211455D0*PPPCORD(III,2,1,3)
               PPPCORD(III,3,1,5) =
     &          0.4436538d0*PPPCORD(III,3,1,1)+0.2352006D0*PPPCORD(III+1,3,1,1)+0.3211455D0*PPPCORD(III,3,1,3)
            ENDDO

            DO III = 1, NMRES
               RES_TYPE = AMINOA(IRES(III))

               IF (III .NE. 1) THEN
                  ATOM_TYPE='N '
                  ID = 4
                  WRITE(27,61)III,ATOM_TYPE,RES_TYPE,III,PPPCORD(III,1,1,ID),PPPCORD(III,2,1,ID),PPPCORD(III,3,1,ID),III
61      FORMAT('ATOM',4X,i3,2X,A2,2X,A3,3X,i3,4X,F8.3,F8.3,F8.3,2X,'1.00',2X,'0.00',6X,'TPDB',1x,I3)

               ENDIF

               ATOM_TYPE='CA'
               ID = 1
               WRITE(27,61)III,ATOM_TYPE,RES_TYPE,III,PPPCORD(III,1,1,ID),PPPCORD(III,2,1,ID),PPPCORD(III,3,1,ID),III

               IF (RES_TYPE .NE. 'gly') THEN
                  ATOM_TYPE='CB'
                  ID = 2
                  WRITE(27,61)III,ATOM_TYPE,RES_TYPE,III,PPPCORD(III,1,1,ID),PPPCORD(III,2,1,ID),PPPCORD(III,3,1,ID),III
               ENDIF

               IF (III .NE. NMRES) THEN
                  ATOM_TYPE='C '
                  ID = 5
                  WRITE(27,61)III,ATOM_TYPE,RES_TYPE,III,PPPCORD(III,1,1,ID),PPPCORD(III,2,1,ID),PPPCORD(III,3,1,ID),III
               ENDIF

               ATOM_TYPE='O '
               ID = 3
               WRITE(27,61)III,ATOM_TYPE,RES_TYPE,III,PPPCORD(III,1,1,ID),PPPCORD(III,2,1,ID),PPPCORD(III,3,1,ID),III

            ENDDO

         ELSE IF (ARNO) THEN
            WRITE(MYUNIT2,'(A,F20.10)') 'N 0.0 0.0 ', 0.577D0
            WRITE(MYUNIT2,'(A,F20.10)') 'O 0.0 0.0 ',-0.577D0
            WRITE(MYUNIT2,65) (QMINP(J1,J2),J2=1,3*(NATOMS-NS))
65          FORMAT('AR ',3F20.10)
         ELSE IF (TOSI.OR.WELCH) THEN
            DO J2=1,NATOMS
               IF (ZSYM(J2).EQ.'PL') WRITE(MYUNIT2,'(A,3F20.10)') 'Na  ',(QMINP(J1,3*(J2-1)+J3),J3=1,3)
               IF (ZSYM(J2).EQ.'MI') WRITE(MYUNIT2,'(A,3F20.10)') 'Cl  ',(QMINP(J1,3*(J2-1)+J3),J3=1,3)
            ENDDO
         ELSE IF (BLJCLUSTER) THEN
            DO J2=1,NATOMS
               IF (J2.LE.NTYPEA) THEN 
                  WRITE(MYUNIT2,'(A,3F20.10)') 'LA  ',(QMINP(J1,3*(J2-1)+J3),J3=1,3)
               ELSE
                  WRITE(MYUNIT2,'(A,3F20.10)') 'LB  ',(QMINP(J1,3*(J2-1)+J3),J3=1,3)
               ENDIF
            ENDDO

         ELSE IF (AMBER) THEN
            DO J2=1,NATOMS
               WRITE(MYUNIT2,'(A,3F20.10)') typech(J2)(1:1),(QMINP(J1,3*(J2-1)+J3),J3=1,3)
            ENDDO
         ELSE IF (AMBERT) THEN
            IF (DUMPSTRUCTURES) THEN
              WRITE(J1CHAR2,'(I3)') J1
              WRITE(J1CHAR,'(A,A)') 'coords.',TRIM(ADJUSTL(J1CHAR2))
              OPEN(UNIT=226,FILE=trim(adjustl(J1CHAR)),STATUS='UNKNOWN')

              DO J2=1,NATOMS
                 WRITE(226,'(3F28.20)') QMINP(J1,3*(J2-1)+1),QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3)
              ENDDO
               WRITE(J1CHAR2,'(I3)') J1
               WRITE(J1CHAR,'(A,A)') 'coords.',TRIM(ADJUSTL(J1CHAR2))
               OPEN(UNIT=226,FILE=trim(adjustl(J1CHAR)),STATUS='UNKNOWN')
              
               DO J2=1,NATOMS
                  WRITE(226,'(3F28.20)') QMINP(J1,3*(J2-1)+1),QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3)
               ENDDO

            ELSE
               DO I1=1,NATOMS
                  WRITE(MYUNIT2,'(A2,3F20.10)') ih(m04+I1-1),COORDS1(3*I1-2),COORDS1(3*I1-1),COORDS1(3*I1)
               ENDDO
            ENDIF

         ELSE IF (ELLIPSOIDT.OR.LJCAPSIDT.OR.GBT.OR.GBDT) THEN
              SUMX=0.0D0
              SUMY=0.0d0
              SUMZ=0.0d0

              DO J2=1,NATOMS/2
                 SUMX=SUMX+QMINP(J1,3*(J2-1)+1)
                 SUMY=SUMY+QMINP(J1,3*(J2-1)+2)
                 SUMZ=SUMZ+QMINP(J1,3*(J2-1)+3)
              ENDDO
              SUMX=2*SUMX/NATOMS
              SUMY=2*SUMY/NATOMS
              SUMZ=2*SUMZ/NATOMS
              DO J2=1,NATOMS/2
                 QMINP(J1,3*(J2-1)+1)=QMINP(J1,3*(J2-1)+1)-SUMX
                 QMINP(J1,3*(J2-1)+2)=QMINP(J1,3*(J2-1)+2)-SUMY
                 QMINP(J1,3*(J2-1)+3)=QMINP(J1,3*(J2-1)+3)-SUMZ
              ENDDO

              DO j2=1,NATOMS/2
                 IF (PARAMONOVPBCX) THEN
                    QMINP(J1,3*J2-2)=QMINP(J1,3*J2-2)-BOXLX*NINT(QMINP(J1,3*J2-2)/BOXLX)
                 ENDIF
                 IF (PARAMONOVPBCY) THEN
                    QMINP(J1,3*J2-1)=QMINP(J1,3*J2-1)-BOXLY*NINT(QMINP(J1,3*J2-1)/BOXLY)
                 ENDIF
                 IF (PARAMONOVPBCZ) THEN
                    QMINP(J1,3*J2  )=QMINP(J1,3*J2  )-BOXLZ*NINT(QMINP(J1,3*J2  )/BOXLZ)
                 ENDIF
              ENDDO

              DO J2=1,NATOMS/2
                 WRITE(MYUNIT2,'(a5,2x,3f20.10,2x,a11,3f20.10)') 'H',QMINP(J1,3*(J2-1)+1), 
     &             QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3),
     &             'atom_vector',QMINP(J1,3*NATOMS/2+3*(J2-1)+1),
     &             QMINP(J1,3*NATOMS/2+3*(J2-1)+2),QMINP(J1,3*NATOMS/2+3*(J2-1)+3)
             ENDDO

         ELSE IF (CHRMMT) THEN
            DO J2=1,NATOMS
               WRITE(MYUNIT2,'(A,1X,3F20.10)') ZSYM(J2)(1:1),(QMINP(J1,3*(J2-1)+J3),J3=1,3)
            ENDDO
            DO J2=1,NATOMS
               DCOORDS(3*(J2-1)+1)=QMINP(J1,3*(J2-1)+1)
               DCOORDS(3*(J2-1)+2)=QMINP(J1,3*(J2-1)+2)
               DCOORDS(3*(J2-1)+3)=QMINP(J1,3*(J2-1)+3)
            ENDDO


         ELSE IF (DBPT .OR. DBPTDT .OR. DMBLMT .OR. LINRODT .OR. LWOTPT .OR. MSTBINT .OR. MSSTOCKT .OR. NCAPT .OR. NPAHT  
     &           .OR. NTIPT .OR. STOCKAAT .OR. PAHAT .OR. PAHW99T .OR. TDHDT .OR. WATERDCT .OR. WATERKZT .OR. PATCHY .OR. PAPT) THEN
            DO J2 = 1, NATOMS/2
               WRITE(MYUNIT2,'(3f20.10)') (QMINP(J1,3*(J2-1)+J3),J3=1,3)
            ENDDO
            DO J2 = 1, NATOMS/2
               WRITE(MYUNIT2,'(3f20.10)') (QMINP(J1,3*NATOMS/2+3*(J2-1)+J3),J3=1,3)
            ENDDO

         ELSE IF (GBT.OR.GBDT.OR.GBDPT.OR.PYGT.OR.PYGDPT.OR.MSGBT.OR.MSPYGT) THEN
            DO J2 = 1, NATOMS/2
               WRITE(MYUNIT2,'(3f20.10)') (QMINP(J1,3*(J2-1)+J3),J3=1,3)    
            ENDDO
            DO J2 = 1, NATOMS/2
               WRITE(MYUNIT2,'(3f20.10)') (QMINP(J1,3*NATOMS/2+3*(J2-1)+J3),J3=1,3)
            ENDDO

         ELSE IF (GEMT) THEN
            DO J2 = 1, NATOMS
               WRITE(MYUNIT2,'(3f20.10)') (QMINP(J1,3*(J2-1)+J3),J3=1,3)
            ENDDO

         ELSE IF (BLNT.AND.(.NOT.P46).AND.(.NOT.G46)) THEN
            DO J2=1,NATOMS
               WRITE(MYUNIT2,'(2A1,1X,3F20.10)') BEADLETTER(J2),'L',(QMINP(J1,3*(J2-1)+J3),J3=1,3)
            ENDDO
         ELSE
            IF (CSMT.AND.(.NOT.SYMMETRIZECSM)) THEN
               WRITE(MYUNIT3,'(I6)') NATOMS
               WRITE(MYUNIT3,'(A,I6,2(A,G20.10))') 'averaged structure for final solution ',J1,
     &                                          ' CSM=',QMINAV(J1),' CSM for reference structure=',QMIN(J1)
               WRITE(MYUNIT3,30) (QMINPCSMAV(J1,J2),J2=1,3*(NATOMS-NS))
            ENDIF
            WRITE(MYUNIT2,30) (QMINP(J1,J2),J2=1,3*(NATOMS-NS))
30          FORMAT('LA ',3F20.10)
         ENDIF

         IF (ASAOOS) THEN
            OPEN(31,file='particles.xyz')
            WRITE(31,*) NATOMS
            WRITE(31,*) ' '
            DO J2=1,NATOMS
                  WRITE(31,'(A,3F20.10)') 'H  ',(QMINP(J1,3*(J2-1)+J3),J3=1,3)
            ENDDO
         END IF

         IF ((NS.GT.0).AND.(.NOT.(WELCH.OR.TOSI))) THEN
            IF (MSORIGT.OR.FRAUSIT) THEN
               WRITE(MYUNIT2,40) (QMINP(J1,J2),J2=3*(NATOMS-NS)+1,3*NATOMS)
40             FORMAT('Si',3F20.10)
            ELSE IF (MSTRANST) THEN
               WRITE(MYUNIT2,40) (QMINP(J1,J2),J2=3*(NATOMS-NS)+1,3*NATOMS)
            ELSE
               WRITE(MYUNIT2,50) (QMINP(J1,J2),J2=3*(NATOMS-NS)+1,3*NATOMS)
50             FORMAT('LB',3F20.10)
            ENDIF
         ENDIF
         IF (AMBER) CALL AMBERDUMP(J1,QMINP)
      ENDDO

      IF (CSMT.AND.(.NOT.SYMMETRIZECSM)) CLOSE(MYUNIT3)
      IF (A9INTET) THEN
         IF (MPIT) THEN
            WRITE (ISTR, '(I10)') MYUNIT-22980+1
            MYUNIT2=(MYUNIT-22980+1)+100
            MYFILENAME2="intelowest."//trim(adjustl(istr))
            OPEN(MYUNIT2,FILE=trim(adjustl(MYFILENAME2)), STATUS="unknown", form="formatted")
         ELSE
            MYUNIT2=25
            OPEN(MYUNIT2,FILE='intelowest',STATUS='UNKNOWN')
         ENDIF
      ENDIF
      IF (A9INTET.AND.AMBERT) THEN
         DO J1=1,NSAVEINTE
            WRITE(MYUNIT2,*) NATOMS
            WRITE(MYUNIT2,10) J1, INTEQMIN(J1), INTEFF(J1)
            IF (DUMPSTRUCTURES) THEN
              WRITE(J1CHAR2,'(I3)') J1
              WRITE(J1CHAR,'(A,A)') 'intecoords.',TRIM(ADJUSTL(J1CHAR2))
              OPEN(UNIT=226,FILE=trim(adjustl(J1CHAR)),STATUS='UNKNOWN')

              DO J2=1,NATOMS
                 WRITE(226,'(3F28.20)') INTEQMINP(J1,3*(J2-1)+1),INTEQMINP(J1,3*(J2-1)+2),INTEQMINP(J1,3*(J2-1)+3)
              ENDDO
               WRITE(J1CHAR2,'(I3)') J1
               WRITE(J1CHAR,'(A,A)') 'intecoords.',TRIM(ADJUSTL(J1CHAR2))
               OPEN(UNIT=226,FILE=trim(adjustl(J1CHAR)),STATUS='UNKNOWN')
              
               DO J2=1,NATOMS
                  WRITE(226,'(3F28.20)') INTEQMINP(J1,3*(J2-1)+1),INTEQMINP(J1,3*(J2-1)+2),INTEQMINP(J1,3*(J2-1)+3)
               ENDDO

            ELSE
               DO I1=1,NATOMS
                  WRITE(MYUNIT2,'(A2,3F20.10)') ih(m04+I1-1),COORDS1(3*I1-2),COORDS1(3*I1-1),COORDS1(3*I1)
               ENDDO
            ENDIF
         ENDDO
      ENDIF

      IF (CHRMMT.AND.RMST) THEN
         OPEN(UNIT=MYUNIT2,FILE='rmsbest',STATUS='UNKNOWN')
         DO J2=1,RMSSAVE
            WRITE(MYUNIT2,'(I6,F6.3,F15.5)')J2,RMSBEST(J2,1),RMSBEST(J2,2)
            WRITE(CRMS,'(I6)') J2
            DCOORDS(1:3*NATOMS)=RMSCOOR(J2,1:3*NATOMS)
            IF(RMSBEST(J2,2)<0.D0) CALL CHARMMDUMP(DCOORDS,'rms.'//TRIM(ADJUSTL(CRMS)))
         ENDDO
      ENDIF

      IF (LJCOULT) THEN
         OPEN(UNIT=26,FILE='ljcoul.xyz',STATUS='UNKNOWN')
         DO J1=1,NSAVE
            WRITE(26,'(I6)') NATOMS
            WRITE(26,10) J1, QMIN(J1), FF(J1)
            DO J2=1,NATOMS
               IF (J2.LE.COULN) THEN
                  WRITE(26,'(A4,3F18.10,A12,3F18.10)') 'O ',
     &                    QMINP(J1,3*(J2-1)+1),QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3)
               ELSE
                  WRITE(26,'(A4,3F18.10,A12,3F18.10)') 'N ',
     &                    QMINP(J1,3*(J2-1)+1),QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3)
               END IF
            ENDDO
         ENDDO


      ELSE IF (STOCKT) THEN
         OPEN(UNIT=26,FILE='stock.xyz',STATUS='UNKNOWN')
         DO J1=1,NSAVE
            WRITE(26,'(I6)') NATOMS/2
            WRITE(26,10) J1, QMIN(J1), FF(J1)
            DO J2=1,NATOMS/2
               WRITE (26,'(A4,3F18.10,A12,3F18.10)') 'LA ',
     &            QMINP(J1,3*(J2-1)+1),QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3),
     &            'atom_vector',
     &            SIN(QMINP(J1,3*(NATOMS/2)+3*(J2-1)+1))*COS(QMINP(J1,3*(NATOMS/2)+3*(J2-1)+2)),
     &            SIN(QMINP(J1,3*(NATOMS/2)+3*(J2-1)+1))*SIN(QMINP(J1,3*(NATOMS/2)+3*(J2-1)+2)),
     &                    COS(QMINP(J1,3*(NATOMS/2)+3*(J2-1)+1))
            ENDDO
         ENDDO

      ELSE IF (ELLIPSOIDT.OR.LJCAPSIDT.OR.GBT.OR.GBDT.OR.MULTISITEPYT) THEN
         DO J1=1,NSAVE
            WRITE(J1CHAR2,'(I3)') J1
          IF (MPIT) THEN
              WRITE (ISTR, '(I10)') MYUNIT-22980+1
              MYUNIT2=(MYUNIT-22980+1)+100
              MYFILENAME2='coords.'//TRIM(ADJUSTL(J1CHAR2))//'.'//trim(adjustl(istr))
              OPEN(MYUNIT2,FILE=trim(adjustl(MYFILENAME2)), STATUS="unknown", form="formatted")
          ELSE 
              MYUNIT2=226
              MYFILENAME2='coords.'//TRIM(ADJUSTL(J1CHAR2))
              OPEN(MYUNIT2,FILE=trim(adjustl(MYFILENAME2)), STATUS="unknown", form="formatted")
          END IF

            DO J2=1,NATOMS
               WRITE(MYUNIT2,'(3F28.20)') QMINP(J1,3*(J2-1)+1),QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3)
            ENDDO


         ENDDO
      IF (MPIT) THEN
         WRITE (ISTR, '(I10)') MYUNIT-22980+1
         MYUNIT2=(MYUNIT-22980+1)+300
         MYFILENAME2="ellipsoid."//trim(adjustl(istr))//".xyz"
         OPEN(MYUNIT2,FILE=trim(adjustl(MYFILENAME2)), STATUS="unknown", form="formatted")
      ELSE
         MYUNIT2=26
         OPEN(MYUNIT2,FILE='ellipsoid.xyz',STATUS='UNKNOWN')
      ENDIF

        do J1=1,NSAVE
        IF(.NOT.MULTISITEPYT) THEN
         WRITE(MYUNIT2,*) NATOMS/2
        ELSE
         WRITE(MYUNIT2,*) NATOMS*NPYSITE/2
        END IF
         WRITE(MYUNIT2,10) J1, QMIN(J1), FF(J1)

            DO J2=1,NATOMS/2
               
               IF (GAYBERNET) THEN
     &                    QMINP(J1,3*NATOMS/2+3*(J2-1)+3),EulerPhi,EulerPsi,EulerTheta,EulerPhiDeg,EulerPsiDeg,EulerThetaDeg)
                  WRITE(MYUNIT2,'(a5,2x,3f20.10,2x,a8,6f15.8,2x,a11,3f15.8)') '0',QMINP(J1,3*(J2-1)+1),
     &                  QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3),
     &                  'ellipse ',1.0D0,1.0D0,GBANISOTROPYR,
     &                                                 EulerPsiDeg,EulerPhiDeg,0.0D0, ! this is in degrees
     &                  'atom_vector',QMINP(J1,3*NATOMS/2+3*(J2-1)+1),
     &                  QMINP(J1,3*NATOMS/2+3*(J2-1)+2),QMINP(J1,3*NATOMS/2+3*(J2-1)+3)
               ELSE IF (PARAMONOVT) THEN
     &                    QMINP(J1,3*NATOMS/2+3*(J2-1)+3),EulerPhi,EulerPsi,EulerTheta,EulerPhiDeg,EulerPsiDeg,EulerThetaDeg)
                  WRITE(MYUNIT2,'(a5,2x,3f20.10,2x,a8,6f15.8,2x,a11,3f15.8)') 'O',QMINP(J1,3*(J2-1)+1),
     &             QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3),
     &             'ellipse ',PARAMc1*2.0D0,PARAMb1*2.0D0,PARAMa1*2.0D0,EulerPhiDeg,EulerPsiDeg,EulerThetaDeg,
     &             'atom_vector',QMINP(J1,3*NATOMS/2+3*(J2-1)+1),
     &             QMINP(J1,3*NATOMS/2+3*(J2-1)+2),QMINP(J1,3*NATOMS/2+3*(J2-1)+3)
               ELSE IF (PYGPERIODICT.OR.PYBINARYT) THEN
     &                           QMINP(J1,3*NATOMS/2+3*(J2-1)+3),EulerPhiDeg,EulerPsiDeg,EulerThetaDeg)

                  WRITE(MYUNIT2,'(a5,2x,3f20.10,2x,a8,6f15.8,2x,a11,3f15.8)') 'O',QMINP(J1,3*(J2-1)+1),
     &                  QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3),
     &            'ellipse ',PYA1BIN(J2,1)*2.0D0,PYA1BIN(J2,2)*2.0D0,PYA1BIN(J2,3)*2.0D0,EulerPhiDeg,EulerPsiDeg,EulerThetaDeg,
     &            'atom_vector',QMINP(J1,3*NATOMS/2+3*(J2-1)+1),
     &            QMINP(J1,3*NATOMS/2+3*(J2-1)+2),QMINP(J1,3*NATOMS/2+3*(J2-1)+3)
               ELSE IF (MULTISITEPYT) THEN

                DO J3=1,NPYSITE

                  WRITE(MYUNIT2,'(a5,2x,3f20.10,2x,a8,12f15.8,2x,a11,3f15.8)') 'O',
     &                                          SITECOORDS(J3,1),SITECOORDS(J3,2),SITECOORDS(J3,3),
     &            'ellipse ',ELLST1(J3,1)*2.0D0,ELLST1(J3,2)*2.0D0,ELLST1(J3,3)*2.0D0,
     &                       ELLMAT(J3,1,1),ELLMAT(J3,1,2),ELLMAT(J3,1,3),ELLMAT(J3,2,1),ELLMAT(J3,2,2),ELLMAT(J3,2,3),
     &                       ELLMAT(J3,3,1),ELLMAT(J3,3,2),ELLMAT(J3,3,3),
     &            'atom_vector',QMINP(J1,3*NATOMS/2+3*(J2-1)+1),
     &            QMINP(J1,3*NATOMS/2+3*(J2-1)+2),QMINP(J1,3*NATOMS/2+3*(J2-1)+3)
                END DO
               ELSE IF (GBT.OR.GBDT) THEN
     &                           QMINP(J1,3*NATOMS/2+3*(J2-1)+3),EulerPhiDeg,EulerPsiDeg,EulerThetaDeg)

                  WRITE(MYUNIT2,'(a5,2x,3f20.10,2x,a8,6f15.8,2x,a11,3f15.8)') 'O',QMINP(J1,3*(J2-1)+1),
     &                  QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3),
     &            'ellipse ',GBKAPPA,1.0D0,1.0D0,EulerPhiDeg,EulerPsiDeg,EulerThetaDeg,
     &            'atom_vector',QMINP(J1,3*NATOMS/2+3*(J2-1)+1),
     &            QMINP(J1,3*NATOMS/2+3*(J2-1)+2),QMINP(J1,3*NATOMS/2+3*(J2-1)+3)
               ELSE IF (LJCAPSIDT) THEN
     &                           QMINP(J1,3*NATOMS/2+3*(J2-1)+3),EulerPhiDeg,EulerPsiDeg,EulerThetaDeg)

                  WRITE(MYUNIT2,'(a5,2x,3f20.10,2x,a8,6f15.8,2x,a11,3f15.8)') 'O',QMINP(J1,3*(J2-1)+1),
     &                  QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3),
     &            'ellipse ',1.0D0,1.0D0,1.0D0,EulerPhiDeg,EulerPsiDeg,EulerThetaDeg,
     &            'atom_vector',QMINP(J1,3*NATOMS/2+3*(J2-1)+1),
     &            QMINP(J1,3*NATOMS/2+3*(J2-1)+2),QMINP(J1,3*NATOMS/2+3*(J2-1)+3)
               ENDIF
            ENDDO
         ENDDO

      ELSE IF (TIP) THEN
         OPEN(UNIT=26,FILE='tip.xyz',STATUS='UNKNOWN')
         DO J1=1,NSAVE
            WRITE(26,'(I6)') (NATOMS/2)*3
            WRITE(26,10) J1, QMIN(J1), FF(J1)
            DO J2=1,NATOMS/2
     1            QMINP(J1,3*(NATOMS/2+J2-1)+1),QMINP(J1,3*(NATOMS/2+J2-1)+2),QMINP(J1,3*(NATOMS/2+J2-1)+3),
     2                     RBCOORDS)
               WRITE(26,'(A4,3F20.10)') 'O ',RBCOORDS(1),RBCOORDS(2),RBCOORDS(3)
               WRITE(26,'(A4,3F20.10)') 'H ',RBCOORDS(4),RBCOORDS(5),RBCOORDS(6)
               WRITE(26,'(A4,3F20.10)') 'H ',RBCOORDS(7),RBCOORDS(8),RBCOORDS(9)
            ENDDO
         ENDDO
      ELSE IF (CAPSID) THEN
         OPEN(UNIT=26,FILE='capsid.xyz',STATUS='UNKNOWN')
         DO J1=1,NSAVE
            WRITE(26,'(I6)') (NATOMS/2)*6
            WRITE(26,10) J1, QMIN(J1), FF(J1)
            DO J2=1,NATOMS/2
     1                       QMINP(J1,3*(NATOMS/2+J2-1)+1),QMINP(J1,3*(NATOMS/2+J2-1)+2),QMINP(J1,3*(NATOMS/2+J2-1)+3),
     2                       RBCOORDS,RAD,HEIGHT)
               DO J3=1,5
                  WRITE(26,'(A4,3F20.10)') 'C1 ',RBCOORDS(3*(J3-1)+1),RBCOORDS(3*(J3-1)+2),RBCOORDS(3*(J3-1)+3)
               ENDDO
               WRITE(26,'(A4,3F20.10)') 'C4  ',RBCOORDS(16),RBCOORDS(17),RBCOORDS(18)
            ENDDO
         ENDDO

      ELSE IF (DBPT) THEN

         OPEN(UNIT=26, FILE='dbp.xyz', STATUS='UNKNOWN')
         GTEST = .FALSE.
         DU    = (/0.D0, 1.D0, 0.D0/) 
         DO J1 = 1, NSAVE

            WRITE(26,'(I6)') (NATOMS/2)*NRBSITES 
            WRITE(26,10) J1, QMIN(J1), FF(J1)

            DO J2 = 1, NATOMS/2

               J3   = 3*J2
               J5   = 3*NATOMS/2 + J3
               P(:) = QMINP(J1,J5-2:J5)


               DO J4 = 1, NRBSITES 

                  IF (J4 == 1) THEN
                     RBCOORDS(1:3) = QMINP(J1,J3-2:J3) + MATMUL(RMI(:,:),SITE(J4,:))
                     WRITE(26,'(A4,3F20.10)') 'O', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                  ELSEIF (J4 == 2) THEN
                     RBCOORDS(1:3) = QMINP(J1,J3-2:J3) + MATMUL(RMI(:,:),SITE(J4,:))
                     WRITE(26,'(A4,3F20.10)') 'C', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                  ELSE
                     RBCOORDS(1:3) = MATMUL(RMI(:,:),DU)
                  WRITE(26,'(A4,3F20.10,2X,A12,2X,3F20.10)')  
     &             'H', QMINP(J1,J3-2), QMINP(J1,J3-1), QMINP(J1,J3),  
     &             'atom_vector', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3) 
                  ENDIF

               ENDDO

            ENDDO

         ENDDO

         RETURN

      ELSE IF (DBPTDT) THEN

         RETURN

      ELSE IF (DMBLMT) THEN

         RETURN

      ELSE IF (GBT .OR. GBDT) THEN

         OPEN(UNIT=26, FILE='gbe.xyz', STATUS='UNKNOWN')
         GTEST = .FALSE.
         DO J1 = 1, NSAVE

            WRITE(26,'(I6)') NATOMS/2
            WRITE(26,10) J1, QMIN(J1), FF(J1)

            DO J2 = 1, NATOMS/2

               J3   = 3*J2
               J5   = 3*NATOMS/2 + J3
               RBCOORDS(1:3) = QMINP(J1,J3-2:J3)
               P(:)          = QMINP(J1,J5-2:J5)


               PHI   = DATAN2(RMI(2,3),RMI(1,3))
               IF (PHI <= 0.D0) PHI = PHI + 2.D0*PI

               THT   = DACOS(RMI(3,3))

               PHI   = PHI*180.D0/PI
               THT   = THT*180.D0/PI

               WRITE(26,'(a5,2x,3f20.10,2x,a8,6f20.10)')
     &         'O', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3),
     &         'ellipse', 2.D0*ESA(1), 2.D0*ESA(2), 2.D0*ESA(3), PHI, THT, 0.D0

            ENDDO

         ENDDO

         RETURN

      ELSE IF (LINRODT) THEN

         OPEN(UNIT=26, FILE='linrod.xyz', STATUS='UNKNOWN')
         GTEST = .FALSE.

         DO J1 = 1, NSAVE

            WRITE(26,'(I6)') (NATOMS/2)*NRBSITES
            WRITE(26,10) J1, QMIN(J1), FF(J1)

            DO J2 = 1, NATOMS/2

               J3    = 3*J2
               J5    = 3*NATOMS/2 + J3
               P(:)  = QMINP(J1,J5-2:J5)



               DO J4 = 1, NRBSITES

                  RBCOORDS(3*J4-2:3*J4) = QMINP(J1,J3-2:J3) + MATMUL(RMI(:,:),SITE(J4,:))

               ENDDO

               DO J4 = 1, NRBSITES

                  J3 = J4 + 1
                  IF (J4 == NRBSITES) J3 = 1
                  P(:) = RBCOORDS(3*J3-2:3*J3) - RBCOORDS(3*J4-2:3*J4)
                  WRITE(26,'(A4,3F20.10,2X,A12,2X,3F20.10)')
     &            'O', RBCOORDS(3*J4-2), RBCOORDS(3*J4-1), RBCOORDS(3*J4), 'atom_vector', P(1), P(2), P(3)

               ENDDO

            ENDDO

         ENDDO


         RETURN

      ELSE IF (LWOTPT) THEN

         OPEN(UNIT=26, FILE='lwotp.xyz', STATUS='UNKNOWN')
         GTEST = .FALSE.
         
         DO J1 = 1, NSAVE

            WRITE(26,'(I6)') (NATOMS/2)*NRBSITES
            WRITE(26,10) J1, QMIN(J1), FF(J1)

            DO J2 = 1, NATOMS/2

               J3    = 3*J2
               J5    = 3*NATOMS/2 + J3
               P(:)  = QMINP(J1,J5-2:J5)



               DO J4 = 1, NRBSITES

                  RBCOORDS(3*J4-2:3*J4) = QMINP(J1,J3-2:J3) + MATMUL(RMI(:,:),SITE(J4,:))

               ENDDO

               DO J4 = 1, NRBSITES

                  J3 = J4 + 1
                  IF (J4 == NRBSITES) J3 = 1
                  P(:) = RBCOORDS(3*J3-2:3*J3) - RBCOORDS(3*J4-2:3*J4)
                  WRITE(26,'(A4,3F20.10,2X,A12,2X,3F20.10)')                 
     &            'O', RBCOORDS(3*J4-2), RBCOORDS(3*J4-1), RBCOORDS(3*J4), 'atom_vector', P(1), P(2), P(3)

               ENDDO

            ENDDO

         ENDDO


         RETURN

      ELSE IF (MSTBINT) THEN

         OPEN(UNIT=26, FILE='mstbin.xyz', STATUS='UNKNOWN')
         GTEST = .FALSE.

         DO J1 = 1, NSAVE

            WRITE(26,'(I6)') NPS*NRBSITES1 + (NATOMS/2 - NPS)*(NRBSITES - NRBSITES1)
            WRITE(26,10) J1, QMIN(J1), FF(J1)

            DO J2 = 1, NATOMS/2

               J3    = 3*J2
               J5    = 3*NATOMS/2 + J3
               P(:)  = QMINP(J1,J5-2:J5)


               IF (J2 <= NPS) THEN
                  NRBS1 = 1
                  NRBS2 = NRBSITES1
               ELSE
                  NRBS1 = NRBSITES1 + 1
                  NRBS2 = NRBSITES
               ENDIF

               DO J4 = NRBS1, NRBS2

                  RBCOORDS(3*J4-2:3*J4) = QMINP(J1,J3-2:J3) + MATMUL(RMI(:,:),SITE(J4,:))

               ENDDO

               DO J4 = NRBS1, NRBS2

                  J3 = J4 + 1
                 
                  IF (J2 <= NPS) THEN
                     IF (J4 == NRBSITES1) J3 = 1
                  ELSE
                     IF (J4 == NRBSITES) J3 = NRBSITES1 + 1
                  ENDIF
 
                  P(:) = RBCOORDS(3*J3-2:3*J3) - RBCOORDS(3*J4-2:3*J4)
                  IF (J2 <= NPS) THEN
                     WRITE(26,'(A4,3F20.10,2X,A12,2X,3F20.10)')
     &               'O', RBCOORDS(3*J4-2), RBCOORDS(3*J4-1), RBCOORDS(3*J4), 'atom_vector', P(1), P(2), P(3)
                  ELSE
                     WRITE(26,'(A4,3F20.10,2X,A12,2X,3F20.10)')
     &               'N', RBCOORDS(3*J4-2), RBCOORDS(3*J4-1), RBCOORDS(3*J4), 'atom_vector', P(1), P(2), P(3)
                  ENDIF

               ENDDO

            ENDDO

         ENDDO


         RETURN

      ELSE IF (MSSTOCKT) THEN

         OPEN(UNIT=26, FILE='msstock.xyz', STATUS='UNKNOWN')
         GTEST = .FALSE.
         
         DO J1 = 1, NSAVE

            WRITE(26,'(I6)') (NATOMS/2)*NRBSITES
            WRITE(26,10) J1, QMIN(J1), FF(J1)

            DO J2 = 1, NATOMS/2

               J3   = 3*J2
               J5   = 3*NATOMS/2 + J3
               P(:) = QMINP(J1,J5-2:J5)


               DO J4 = 1, NRBSITES

                  RBCOORDS(1:3) = QMINP(J1,J3-2:J3) + MATMUL(RMI(:,:),SITE(J4,:))
                  P(:)          = MATMUL(RMI(:,:),RBUV(J4,:))
                  WRITE(26,'(A4,3F20.10,2X,A12,2X,3F20.10)')                          
     &            'O', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3), 'atom_vector', P(1), P(2), P(3)

               ENDDO

            ENDDO

         ENDDO


         OPEN(UNIT=28, FILE='msstktr.xyz', STATUS='UNKNOWN')
         GTEST = .FALSE.

         DO J1 = 1, NSAVE

            WRITE(28,'(I6)') (NATOMS/2)*NRBSITES !(NRBSITES - 1)
            WRITE(28,10) J1, QMIN(J1), FF(J1)

            DO J2 = 1, NATOMS/2

               J3    = 3*J2
               J5    = 3*NATOMS/2 + J3
               P(:)  = QMINP(J1,J5-2:J5)


               DO J4 = 1, NRBSITES

                  RBCOORDS(3*J4-2:3*J4) = QMINP(J1,J3-2:J3) + MATMUL(RMI(:,:),SITE(J4,:))

               ENDDO

               DO J4 = 1, NRBSITES !- 1

                  J3 = J4 + 1
                  IF (J4 == NRBSITES) J3 = 1
                  P(:) = RBCOORDS(3*J3-2:3*J3) - RBCOORDS(3*J4-2:3*J4)
                  WRITE(28,'(A4,3F20.10,2X,A12,2X,3F20.10)')                 
     &            'O', RBCOORDS(3*J4-2), RBCOORDS(3*J4-1), RBCOORDS(3*J4), 'atom_vector', P(1), P(2), P(3)

               ENDDO

            ENDDO

         ENDDO


         RETURN

      ELSE IF (MULTPAHAT) THEN

         RETURN

      ELSE IF (NPAHT .OR. PAHAT .OR. PAHW99T) THEN

         OPEN(UNIT=26, FILE='rigid.xyz', STATUS='UNKNOWN')
         GTEST = .FALSE.

         IF (PAHW99T) THEN
            NCPHST = NCARBON + (NRBSITES-NCARBON)/2
            DO J1 = 1, (NCPHST-NCARBON) 
               SITE(NCARBON+J1,:) = SITE(NCPHST+J1,:)
            ENDDO
         ELSE
            NCPHST = NRBSITES
         ENDIF

         DO J1 = 1, NSAVE

            WRITE(26,'(I6)') (NATOMS/2)*NRBSITES
            WRITE(26,10) J1, QMIN(J1), FF(J1)

            DO J2 = 1, NATOMS/2

               J3   = 3*J2
               J5   = 3*NATOMS/2 + J3
               P(:) = QMINP(J1,J5-2:J5)


               DO J4 = 1, NCPHST

                  RBCOORDS(1:3) = QMINP(J1,J3-2:J3) + MATMUL(RMI(:,:),SITE(J4,:))
                  IF (J4 <= NCARBON) THEN
                     WRITE(26,'(A4,3F20.10)') 'C', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                  ELSE
                     WRITE(26,'(A4,3F20.10)') 'H', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                  ENDIF

               ENDDO

            ENDDO

         ENDDO

         RETURN

      ELSE IF (NCAPT) THEN

         RETURN

      ELSE IF (NTIPT) THEN

         RETURN

      ELSE IF (PAPT) THEN

         RETURN


      ELSE IF (PATCHY) THEN

         RETURN


      ELSE IF (PYGT .OR. PYGDPT) THEN

         OPEN(UNIT=26, FILE='ellipsoids.xyz', STATUS='UNKNOWN')
         GTEST = .FALSE.

         DO J1 = 1, NSAVE

            WRITE(26,'(I6)') (NATOMS/2)
            WRITE(26,10) J1, QMIN(J1), FF(J1)

            DO J2 = 1, NATOMS/2

               J3   = 3*J2
               J5   = 3*NATOMS/2 + J3

               RBCOORDS(1:3) = QMINP(J1,J3-2:J3)
               P(:)          = QMINP(J1,J5-2:J5)


               PHI = DATAN2(RMI(2,3),RMI(1,3))
               IF (PHI <= 0.D0) PHI = PHI + 2.D0*PI

               THT = DACOS(RMI(3,3))

               IF (CHI <= 0.D0) CHI = CHI + 2.D0*PI

               PHI = PHI*180.D0/PI
               THT = THT*180.D0/PI

               WRITE(26,'(a5,2x,3f20.10,2x,a8,6f20.10)')                 
     &              'O', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3),         
     &              'ellipse', 2.D0*PYA1(1), 2.D0*PYA1(2), 2.D0*PYA1(3), PHI, THT, CHI 

            ENDDO

         ENDDO


         RETURN

      ELSE IF (SILANET) THEN

         RETURN

      ELSE IF (STOCKAAT) THEN

         OPEN(UNIT=26, FILE='stockaa.xyz', STATUS='UNKNOWN')
         GTEST = .FALSE.
         DU    = (/0.D0, 0.D0, 1.D0/)

         DO J1 = 1, NSAVE

            WRITE(26,'(I6)') NATOMS/2
            WRITE(26,10) J1, QMIN(J1), FF(J1)

            DO J2 = 1, NATOMS/2

               J3   = 3*J2
               J5   = 3*NATOMS/2 + J3
               P(:) = QMINP(J1,J5-2:J5)


               RBCOORDS(1:3) = MATMUL(RMI(:,:),DU(:))
               WRITE(26,'(A4,3F20.10,2X,A12,2X,3F20.10)') 'O', QMINP(J1,J3-2), QMINP(J1,J3-1), QMINP(J1,J3), 
     &         'atom_vector', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)


            ENDDO

         ENDDO


         RETURN

      ELSE IF (TDHDT) THEN

         RETURN

      ELSE IF (WATERDCT .OR. WATERKZT) THEN

         OPEN(UNIT=26, FILE='rigid.xyz', STATUS='UNKNOWN')
         GTEST = .FALSE.

         DO J1 = 1, NSAVE

            WRITE(26,'(I6)') (NATOMS/2)*(NRBSITES - 1)
            WRITE(26,10) J1, QMIN(J1), FF(J1)

            DO J2 = 1, NATOMS/2

               J3   = 3*J2
               J5   = 3*NATOMS/2 + J3
               P(:) = QMINP(J1,J5-2:J5)


               DO J4 = 1, NRBSITES - 1

                  RBCOORDS(1:3) = QMINP(J1,J3-2:J3) + MATMUL(RMI(:,:),SITE(J4,:))
                  IF (J4 == 1) THEN
                     WRITE(26,'(A4,3F20.10)') 'O', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                  ELSE
                     WRITE(26,'(A4,3F20.10)') 'H', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                  ENDIF

               ENDDO

            ENDDO

         ENDDO

         RETURN

      ELSE IF (RIGID) THEN
         OPEN(UNIT=26,FILE='rigid.xyz',STATUS='UNKNOWN')
         DO J1=1,NSAVE
            WRITE(26,'(I6)') (NATOMS/2)*NRBSITES
            WRITE(26,10) J1, QMIN(J1), FF(J1)
            DO J2=1,NATOMS/2
     1                       QMINP(J1,3*(NATOMS/2+J2-1)+1),
     2                       QMINP(J1,3*(NATOMS/2+J2-1)+2),
     3                       QMINP(J1,3*(NATOMS/2+J2-1)+3),
     4                       RBCOORDS,NRBSITES,SITE)
               DO J3=1,NRBSITES
                  WRITE(26,'(A4,3F20.10)') 'LA ',RBCOORDS(3*(J3-1)+1),RBCOORDS(3*(J3-1)+2),RBCOORDS(3*(J3-1)+3)
               ENDDO
            ENDDO
         ENDDO

      ENDIF

      RETURN
      END

      include "finalio.amberdump.i.f"
      include "finalio.capsidio.i.f"
      include "finalio.rbio.i.f"
      include "finalio.tipio.i.f"
