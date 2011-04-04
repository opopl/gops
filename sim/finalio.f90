!op226>=================================== 
!op226> GPL License Info {{{ 
C   GMIN: A program for finding global minima
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of GMIN.
C
C   GMIN is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   GMIN is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
!op226>}}} 
!op226>=================================== 
      SUBROUTINE FINALIO
!op226>=================================== 
!> \brief Produce final quenches  
!op226>=================================== 
!op226> Declarations {{{ 
      USE commons
      use modamber
      use modamber9, only : coords1,lcrd,ih,m04,natom,ambfinalio_node
      use pymodule, only : SITECOORDS,ELLST1,ELLMAT
      use qmodule
      USE modcharmm
      USE AMHGLOBALS, ONLY:NMRES,IRES

      IMPLICIT NONE

C   MCP
      INTEGER III, I3,  GLY_COUNT, ID, NUMCRD, NUMPRO, NCPHST
      INTEGER J1, J2, J3, J4, J5, NTYPEA, MYUNIT2, I1, NDUMMY, MYUNIT3, NC, NRBS1, NRBS2
      DOUBLE PRECISION EPSAB, EPSBB, SIGAB, SIGBB, RBCOORDS(NRBSITES*3), DCOORDS(3*NATOMS)
      DOUBLE PRECISION P3(3,3), P(3), DU(3), RMI(3,3), DRMI(3,3), PI, PHI, THT, CHI
      DOUBLE PRECISION, ALLOCATABLE :: XCOORDS(:), YCOORDS(:)
      COMMON /BIN/ EPSAB, EPSBB, SIGAB, SIGBB, NTYPEA
      CHARACTER(LEN=13) J1CHAR,J1CHAR2                  !for gay-berne output files
      DOUBLE PRECISION EulerPhi,EulerPsi,EulerTheta,EulerThetadeg,EulerPhiDeg,EulerPsiDeg  ! Euler angles for ellipsoids of revolution

      DOUBLE PRECISION EPS2, RAD, HEIGHT,sumx,sumy,sumz
      LOGICAL :: GTEST
      COMMON /CAPS/ EPS2, RAD, HEIGHT

      CHARACTER(LEN=6) :: CRMS
      CHARACTER(LEN=20) :: MYFILENAME2, ISTR, DBNUM, MYFILENAME3
      CHARACTER(LEN=15), ALLOCATABLE :: DBNAME(:) 

C  AMH 
      CHARACTER(LEN=3) :: RES_TYPE
      CHARACTER(LEN=2) :: ATOM_TYPE
      CHARACTER*1 COUNTTT
      INTEGER COUNTT
      DOUBLE PRECISION  PPPCORD(NMRES*3*3,3,3,5)
      EXTERNAL NUM_TO_CHAR
!op226> End declarations }}} 

      PI = 4.D0*DATAN(1.D0)

      NUMPRO = 1
      NUMCRD = 3

      ALLOCATE(DBNAME(NSAVE))
     
      DO J1=1,NSAVE
         WRITE(DBNUM,*) J1
         DBNAME(J1)='dbase.'//TRIM(ADJUSTL(DBNUM))
      ENDDO
      DCOORDS(1:3*NATOMS)=0.0D0
!      IF (AMH) THEN
!         CALL WALESAMH_FINALIO
!         STOP
!      ENDIF

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
            COUNTT=J1
            CALL NUM_TO_CHAR(COUNTT,COUNTTT)
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
!        IF (CSMT.AND.DEBUG) WRITE(MYUNIT,'(A,I6,2G20.10)') 'finalio> J1,QMIN,QMINAV=',J1,QMIN(J1),QMINAV(J1)
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
C
C   OUTPUT CORDS FOR LOWEST IN X,Y,Z FORMAT
C   OUTPUT CORDS FOR MOVIE_GMIN IN MOVIESEG FORMAT
C
            WRITE(26,683)NUMCRD,j1,NUMCRD,REAL(NUMCRD),NUMCRD
683         FORMAT(3(I6,1X),F8.4,1X,I5,' STUCT SNAP T T TID')

            GLY_COUNT = 0

            DO 1964 III = 1,NMRES
               IF (IRES(III).EQ.8) THEN
!!                pppcord(residue, xyz, numpro, atom types
                  PPPCORD(III, 1, 1, 1) = REAL(QMINP(J1,9*(III-1)+1- GLY_COUNT*3)) !  CA X
                  PPPCORD(III, 2, 1, 1) = REAL(QMINP(J1,9*(III-1)+2- GLY_COUNT*3)) !  CA Y
                  PPPCORD(III, 3, 1, 1) = REAL(QMINP(J1,9*(III-1)+3- GLY_COUNT*3)) !  CA Z
!    SWAP  CA for CB 
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
!632           FORMAT('CA: ',3(f25.15,1x),'CB: ',3(f25.15,1x),'Ox: ', 3(f25.15,1x))
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
            CLOSE(27)

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
          ! sf344> write out coordinates
            COORDS1(1:3*NATOMS) = QMINP(J1,1:3*NATOMS)
            IF (DUMPSTRUCTURES) THEN
              CALL AMBERFINALIO(j1,25,AMBFINALIO_NODE,'0',0,COORDS1(1:3*NATOMS))         
              WRITE(J1CHAR2,'(I3)') J1
              WRITE(J1CHAR,'(A,A)') 'coords.',TRIM(ADJUSTL(J1CHAR2))
              OPEN(UNIT=226,FILE=trim(adjustl(J1CHAR)),STATUS='UNKNOWN')

              DO J2=1,NATOMS
                 WRITE(226,'(3F28.20)') QMINP(J1,3*(J2-1)+1),QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3)
              ENDDO
              CLOSE(226)
               WRITE(J1CHAR2,'(I3)') J1
               WRITE(J1CHAR,'(A,A)') 'coords.',TRIM(ADJUSTL(J1CHAR2))
               OPEN(UNIT=226,FILE=trim(adjustl(J1CHAR)),STATUS='UNKNOWN')
              
               DO J2=1,NATOMS
                  WRITE(226,'(3F28.20)') QMINP(J1,3*(J2-1)+1),QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3)
               ENDDO
               CLOSE(226)

            ELSE
               DO I1=1,NATOMS
                  WRITE(MYUNIT2,'(A2,3F20.10)') ih(m04+I1-1),COORDS1(3*I1-2),COORDS1(3*I1-1),COORDS1(3*I1)
               ENDDO
            ENDIF

         ELSE IF (ELLIPSOIDT.OR.LJCAPSIDT.OR.GBT.OR.GBDT) THEN
C
C determine the centre of coordinates, then centre them
C
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
                        ! ensure y component of particle 1 vector is within BoxLy/2 of zero. 
                        ! If it isn't then subtract integer number of boxly's such that it is.
                    QMINP(J1,3*J2-2)=QMINP(J1,3*J2-2)-BOXLX*NINT(QMINP(J1,3*J2-2)/BOXLX)
                 ENDIF
                 IF (PARAMONOVPBCY) THEN
                        ! ensure y component of particle 1 vector is within BoxLy/2 of zero. 
                        ! If it isn't then subtract integer number of boxly's such that it is.
                    QMINP(J1,3*J2-1)=QMINP(J1,3*J2-1)-BOXLY*NINT(QMINP(J1,3*J2-1)/BOXLY)
                 ENDIF
                 IF (PARAMONOVPBCZ) THEN
                        ! ensure y component of particle 1 vector is within BoxLy/2 of zero. 
                        ! If it isn't then subtract integer number of boxly's such that it is.
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
!       csw34> This DO loop appeared to be be missing on 30/9/08 which would easily
!       explain the output problems! 
            DO J2=1,NATOMS
               DCOORDS(3*(J2-1)+1)=QMINP(J1,3*(J2-1)+1)
               DCOORDS(3*(J2-1)+2)=QMINP(J1,3*(J2-1)+2)
               DCOORDS(3*(J2-1)+3)=QMINP(J1,3*(J2-1)+3)
            ENDDO
            CALL CHARMMDUMP(DCOORDS,DBNAME(J1))

!    DC430 >
!    |gd351> added patchy

         ELSE IF (DBPT .OR. DBPTDT .OR. DMBLMT .OR. LINRODT .OR. LWOTPT .OR. MSTBINT .OR. MSSTOCKT .OR. NCAPT .OR. NPAHT  
     &           .OR. NTIPT .OR. STOCKAAT .OR. PAHAT .OR. PAHW99T .OR. TDHDT .OR. WATERDCT .OR. WATERKZT .OR. PATCHY) THEN
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
C
C this writes 'lowest' in xyz (Xmakemol) format
C
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

!|gd351>
         IF (ASAOOS) THEN
            OPEN(31,file='particles.xyz')
            WRITE(31,*) NATOMS
            WRITE(31,*) ' '
            DO J2=1,NATOMS
                  WRITE(31,'(A,3F20.10)') 'H  ',(QMINP(J1,3*(J2-1)+J3),J3=1,3)
            ENDDO
            CLOSE(31)
         END IF
!<gd351|

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

C
C  End of loop over dump to file lowest or equivalent.
C
      CLOSE(MYUNIT2)
      IF (CSMT.AND.(.NOT.SYMMETRIZECSM)) CLOSE(MYUNIT3)
C
C     csw34> New loop for dumping interaction energy files if A9INTE is specified
C     Added the missing IF block to test for A9INTE 9/12/09 DJW
C
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
C
C     csw34> loop structure copied from the ELSEIF(AMBERT) block above
C
      IF (A9INTET.AND.AMBERT) THEN
         DO J1=1,NSAVEINTE
            WRITE(MYUNIT2,*) NATOMS
!     csw34> write header to intelowest for current minimum
            WRITE(MYUNIT2,10) J1, INTEQMIN(J1), INTEFF(J1)
!     sf344> write out coordinates
            COORDS1(1:3*NATOMS) = INTEQMINP(J1,1:3*NATOMS)
            IF (DUMPSTRUCTURES) THEN
              CALL INTEFINALIO(j1,MYUNIT2,AMBFINALIO_NODE,'0',0,COORDS1(1:3*NATOMS))         
              WRITE(J1CHAR2,'(I3)') J1
              WRITE(J1CHAR,'(A,A)') 'intecoords.',TRIM(ADJUSTL(J1CHAR2))
              OPEN(UNIT=226,FILE=trim(adjustl(J1CHAR)),STATUS='UNKNOWN')

              DO J2=1,NATOMS
                 WRITE(226,'(3F28.20)') INTEQMINP(J1,3*(J2-1)+1),INTEQMINP(J1,3*(J2-1)+2),INTEQMINP(J1,3*(J2-1)+3)
              ENDDO
              CLOSE(226)
               WRITE(J1CHAR2,'(I3)') J1
               WRITE(J1CHAR,'(A,A)') 'intecoords.',TRIM(ADJUSTL(J1CHAR2))
               OPEN(UNIT=226,FILE=trim(adjustl(J1CHAR)),STATUS='UNKNOWN')
              
               DO J2=1,NATOMS
                  WRITE(226,'(3F28.20)') INTEQMINP(J1,3*(J2-1)+1),INTEQMINP(J1,3*(J2-1)+2),INTEQMINP(J1,3*(J2-1)+3)
               ENDDO
               CLOSE(226)

            ELSE
               DO I1=1,NATOMS
                  WRITE(MYUNIT2,'(A2,3F20.10)') ih(m04+I1-1),COORDS1(3*I1-2),COORDS1(3*I1-1),COORDS1(3*I1)
               ENDDO
            ENDIF
         ENDDO
      ENDIF
C
C  End of loop over dump to file intelowest 
C
      CLOSE(MYUNIT2)

!     csw34> Edits to the RMS keyword          
      IF (CHRMMT.AND.RMST) THEN
!        IF (PROGRESS) THEN
!           DCOORDS(1:3*NATOMS)=RMSCOOR(1,1:3*NATOMS)
!           IF(RMSBEST(1,2)<0.D0) CALL CHARMMDUMP(DCOORDS,'closestrms')
!           WRITE(MYUNIT,'(A9,F8.5)') 'RMSDmin= ',RMSBEST(1,1)
!        ELSE      
! REMEMBER TO RE-INDENT THE BELOW IF UNCOMMENTING ABOVE!
         OPEN(UNIT=MYUNIT2,FILE='rmsbest',STATUS='UNKNOWN')
         DO J2=1,RMSSAVE
            WRITE(MYUNIT2,'(I6,F6.3,F15.5)')J2,RMSBEST(J2,1),RMSBEST(J2,2)
            WRITE(CRMS,'(I6)') J2
            DCOORDS(1:3*NATOMS)=RMSCOOR(J2,1:3*NATOMS)
            IF(RMSBEST(J2,2)<0.D0) CALL CHARMMDUMP(DCOORDS,'rms.'//TRIM(ADJUSTL(CRMS)))
         ENDDO
         CLOSE(MYUNIT2)
!        ENDIF
      ENDIF

      IF (LJCOULT) THEN
         OPEN(UNIT=26,FILE='ljcoul.xyz',STATUS='UNKNOWN')
         DO J1=1,NSAVE
            WRITE(26,'(I6)') NATOMS
            WRITE(26,10) J1, QMIN(J1), FF(J1)
            DO J2=1,NATOMS
C              Use "O" atom type to highlight charged particles and "N" for the neutral ones.
               IF (J2.LE.COULN) THEN
                  WRITE(26,'(A4,3F18.10,A12,3F18.10)') 'O ',
     &                    QMINP(J1,3*(J2-1)+1),QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3)
               ELSE
                  WRITE(26,'(A4,3F18.10,A12,3F18.10)') 'N ',
     &                    QMINP(J1,3*(J2-1)+1),QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3)
               END IF
            ENDDO
         ENDDO
         CLOSE(26)


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
         CLOSE(26)

      ELSE IF (ELLIPSOIDT.OR.LJCAPSIDT.OR.GBT.OR.GBDT.OR.MULTISITEPYT) THEN
         DO J1=1,NSAVE
            WRITE(J1CHAR2,'(I3)') J1
!            WRITE(J1CHAR,'(A5,A,A4)') 'gbmin',TRIM(ADJUSTL(J1CHAR2)),'.xyz'
!            OPEN(UNIT=26,FILE=trim(adjustl(J1CHAR)),STATUS='UNKNOWN')
!
!            DO J2=1,NATOMS/2
!               WRITE(26,'(3F18.10)') QMINP(J1,3*(J2-1)+1),QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3)
!            ENDDO
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

            CLOSE(MYUNIT2)

         ENDDO
C
C  Write out lowest NSAVE structures to xmakemol ellipsoid format for 
C  clusters of ellipsoids of revolution 
C
C         OPEN(UNIT=26,FILE="ellipsoid.xyz",STATUS='UNKNOWN')
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
                  CALL EllipsoidsAAtoPolar(QMINP(J1,3*NATOMS/2+3*(J2-1)+1),QMINP(J1,3*NATOMS/2+3*(J2-1)+2),
     &                    QMINP(J1,3*NATOMS/2+3*(J2-1)+3),EulerPhi,EulerPsi,EulerTheta,EulerPhiDeg,EulerPsiDeg,EulerThetaDeg)
!                 EulerPhiDeg = 90-EulerPhiDeg    ! EulerPhiDeg returned from EllipsoidsAAtoPolar is in fact the alpha angle
!                                                ! defined by me (angle of vector with the xy plane)
                  WRITE(MYUNIT2,'(a5,2x,3f20.10,2x,a8,6f15.8,2x,a11,3f15.8)') '0',QMINP(J1,3*(J2-1)+1),
     &                  QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3),
     &                  'ellipse ',1.0D0,1.0D0,GBANISOTROPYR,
     &                                                 EulerPsiDeg,EulerPhiDeg,0.0D0, ! this is in degrees
     &                  'atom_vector',QMINP(J1,3*NATOMS/2+3*(J2-1)+1),
     &                  QMINP(J1,3*NATOMS/2+3*(J2-1)+2),QMINP(J1,3*NATOMS/2+3*(J2-1)+3)
               ELSE IF (PARAMONOVT) THEN
                  CALL EllipsoidsAAtoPolar(QMINP(J1,3*NATOMS/2+3*(J2-1)+1),QMINP(J1,3*NATOMS/2+3*(J2-1)+2),
     &                    QMINP(J1,3*NATOMS/2+3*(J2-1)+3),EulerPhi,EulerPsi,EulerTheta,EulerPhiDeg,EulerPsiDeg,EulerThetaDeg)
                  WRITE(MYUNIT2,'(a5,2x,3f20.10,2x,a8,6f15.8,2x,a11,3f15.8)') 'O',QMINP(J1,3*(J2-1)+1),
     &             QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3),
     &             'ellipse ',PARAMc1*2.0D0,PARAMb1*2.0D0,PARAMa1*2.0D0,EulerPhiDeg,EulerPsiDeg,EulerThetaDeg,
     &             'atom_vector',QMINP(J1,3*NATOMS/2+3*(J2-1)+1),
     &             QMINP(J1,3*NATOMS/2+3*(J2-1)+2),QMINP(J1,3*NATOMS/2+3*(J2-1)+3)
               ELSE IF (PYGPERIODICT.OR.PYBINARYT) THEN
                  CALL AAtoEuler(QMINP(J1,3*NATOMS/2+3*(J2-1)+1),QMINP(J1,3*NATOMS/2+3*(J2-1)+2),
     &                           QMINP(J1,3*NATOMS/2+3*(J2-1)+3),EulerPhiDeg,EulerPsiDeg,EulerThetaDeg)

                  WRITE(MYUNIT2,'(a5,2x,3f20.10,2x,a8,6f15.8,2x,a11,3f15.8)') 'O',QMINP(J1,3*(J2-1)+1),
     &                  QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3),
     &            'ellipse ',PYA1BIN(J2,1)*2.0D0,PYA1BIN(J2,2)*2.0D0,PYA1BIN(J2,3)*2.0D0,EulerPhiDeg,EulerPsiDeg,EulerThetaDeg,
     &            'atom_vector',QMINP(J1,3*NATOMS/2+3*(J2-1)+1),
     &            QMINP(J1,3*NATOMS/2+3*(J2-1)+2),QMINP(J1,3*NATOMS/2+3*(J2-1)+3)
               ELSE IF (MULTISITEPYT) THEN
                 CALL AAtoSites(QMINP(J1,3*(J2-1)+1:3*(J2-1)+3),QMINP(J1,3*NATOMS/2+3*(J2-1)+1:3*NATOMS/2+3*(J2-1)+3),SITECOORDS)

!                  CALL AAtoEuler(QMINPS(J1,3*NPYSITE*NATOMS/2+3*(J2-1)+1),QMINP(J1,3*NPYSITE*NATOMS/2+3*(J2-1)+2),
!     &                           QMINP(J1,3*NPYSITE*NATOMS/2+3*(J2-1)+3),EulerPhiDeg,EulerPsiDeg,EulerThetaDeg)
                DO J3=1,NPYSITE
!                   Call AAtoEuler(ELLST3(J3,1),ELLST3(J3,2),ELLST3(J3,3),EulerPhiDeg,EulerThetaDeg,EulerPsiDeg)

!                  WRITE(MYUNIT2,'(a5,2x,3f20.10)') 'LA',SITECOORDS(J3,1), SITECOORDS(J3,2), SITECOORDS(J3,3)
                  WRITE(MYUNIT2,'(a5,2x,3f20.10,2x,a8,12f15.8,2x,a11,3f15.8)') 'O',
     &                                          SITECOORDS(J3,1),SITECOORDS(J3,2),SITECOORDS(J3,3),
     &            'ellipse ',ELLST1(J3,1)*2.0D0,ELLST1(J3,2)*2.0D0,ELLST1(J3,3)*2.0D0,
     &                       ELLMAT(J3,1,1),ELLMAT(J3,1,2),ELLMAT(J3,1,3),ELLMAT(J3,2,1),ELLMAT(J3,2,2),ELLMAT(J3,2,3),
     &                       ELLMAT(J3,3,1),ELLMAT(J3,3,2),ELLMAT(J3,3,3),
     &            'atom_vector',QMINP(J1,3*NATOMS/2+3*(J2-1)+1),
     &            QMINP(J1,3*NATOMS/2+3*(J2-1)+2),QMINP(J1,3*NATOMS/2+3*(J2-1)+3)
                END DO
               ELSE IF (GBT.OR.GBDT) THEN
                  CALL AAtoEuler(QMINP(J1,3*NATOMS/2+3*(J2-1)+1),QMINP(J1,3*NATOMS/2+3*(J2-1)+2),
     &                           QMINP(J1,3*NATOMS/2+3*(J2-1)+3),EulerPhiDeg,EulerPsiDeg,EulerThetaDeg)

                  WRITE(MYUNIT2,'(a5,2x,3f20.10,2x,a8,6f15.8,2x,a11,3f15.8)') 'O',QMINP(J1,3*(J2-1)+1),
     &                  QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3),
     &            'ellipse ',GBKAPPA,1.0D0,1.0D0,EulerPhiDeg,EulerPsiDeg,EulerThetaDeg,
     &            'atom_vector',QMINP(J1,3*NATOMS/2+3*(J2-1)+1),
     &            QMINP(J1,3*NATOMS/2+3*(J2-1)+2),QMINP(J1,3*NATOMS/2+3*(J2-1)+3)
               ELSE IF (LJCAPSIDT) THEN
                  CALL AAtoEuler(QMINP(J1,3*NATOMS/2+3*(J2-1)+1),QMINP(J1,3*NATOMS/2+3*(J2-1)+2),
     &                           QMINP(J1,3*NATOMS/2+3*(J2-1)+3),EulerPhiDeg,EulerPsiDeg,EulerThetaDeg)

                  WRITE(MYUNIT2,'(a5,2x,3f20.10,2x,a8,6f15.8,2x,a11,3f15.8)') 'O',QMINP(J1,3*(J2-1)+1),
     &                  QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3),
     &            'ellipse ',1.0D0,1.0D0,1.0D0,EulerPhiDeg,EulerPsiDeg,EulerThetaDeg,
     &            'atom_vector',QMINP(J1,3*NATOMS/2+3*(J2-1)+1),
     &            QMINP(J1,3*NATOMS/2+3*(J2-1)+2),QMINP(J1,3*NATOMS/2+3*(J2-1)+3)
               ENDIF
            ENDDO
         ENDDO
         CLOSE(MYUNIT2)

      ELSE IF (TIP) THEN
         OPEN(UNIT=26,FILE='tip.xyz',STATUS='UNKNOWN')
         DO J1=1,NSAVE
            WRITE(26,'(I6)') (NATOMS/2)*3
            WRITE(26,10) J1, QMIN(J1), FF(J1)
            DO J2=1,NATOMS/2
               CALL TIPIO(QMINP(J1,3*(J2-1)+1),QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3),
     1            QMINP(J1,3*(NATOMS/2+J2-1)+1),QMINP(J1,3*(NATOMS/2+J2-1)+2),QMINP(J1,3*(NATOMS/2+J2-1)+3),
     2                     RBCOORDS)
               WRITE(26,'(A4,3F20.10)') 'O ',RBCOORDS(1),RBCOORDS(2),RBCOORDS(3)
               WRITE(26,'(A4,3F20.10)') 'H ',RBCOORDS(4),RBCOORDS(5),RBCOORDS(6)
               WRITE(26,'(A4,3F20.10)') 'H ',RBCOORDS(7),RBCOORDS(8),RBCOORDS(9)
            ENDDO
         ENDDO
         CLOSE(26)
      ELSE IF (CAPSID) THEN
         OPEN(UNIT=26,FILE='capsid.xyz',STATUS='UNKNOWN')
         DO J1=1,NSAVE
            WRITE(26,'(I6)') (NATOMS/2)*6
            WRITE(26,10) J1, QMIN(J1), FF(J1)
            DO J2=1,NATOMS/2
               CALL CAPSIDIO(QMINP(J1,3*(J2-1)+1),QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3),
     1                       QMINP(J1,3*(NATOMS/2+J2-1)+1),QMINP(J1,3*(NATOMS/2+J2-1)+2),QMINP(J1,3*(NATOMS/2+J2-1)+3),
     2                       RBCOORDS,RAD,HEIGHT)
               DO J3=1,5
                  WRITE(26,'(A4,3F20.10)') 'C1 ',RBCOORDS(3*(J3-1)+1),RBCOORDS(3*(J3-1)+2),RBCOORDS(3*(J3-1)+3)
               ENDDO
               WRITE(26,'(A4,3F20.10)') 'C4  ',RBCOORDS(16),RBCOORDS(17),RBCOORDS(18)
            ENDDO
         ENDDO
         CLOSE(26)

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

               CALL RMDRVT(P, RMI, DRMI, DRMI, DRMI, GTEST)

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

         CALL VIEWDMBLTD()
         RETURN

      ELSE IF (DMBLMT) THEN

         CALL VIEWDMBL()
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

               CALL RMDRVT(P, RMI, DRMI, DRMI, DRMI, GTEST)

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

               CALL RMDRVT(P, RMI, DRMI, DRMI, DRMI, GTEST)


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

         CLOSE(UNIT=26)

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

               CALL RMDRVT(P, RMI, DRMI, DRMI, DRMI, GTEST)


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

         CLOSE(UNIT=26)

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

               CALL RMDRVT(P, RMI, DRMI, DRMI, DRMI, GTEST)

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

         CLOSE(UNIT=26)

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

               CALL RMDRVT(P, RMI, DRMI, DRMI, DRMI, GTEST)

               DO J4 = 1, NRBSITES

                  RBCOORDS(1:3) = QMINP(J1,J3-2:J3) + MATMUL(RMI(:,:),SITE(J4,:))
                  P(:)          = MATMUL(RMI(:,:),RBUV(J4,:))
                  WRITE(26,'(A4,3F20.10,2X,A12,2X,3F20.10)')                          
     &            'O', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3), 'atom_vector', P(1), P(2), P(3)

               ENDDO

            ENDDO

         ENDDO

         CLOSE(UNIT=26)

         OPEN(UNIT=28, FILE='msstktr.xyz', STATUS='UNKNOWN')
         GTEST = .FALSE.

         DO J1 = 1, NSAVE

            WRITE(28,'(I6)') (NATOMS/2)*NRBSITES !(NRBSITES - 1)
            WRITE(28,10) J1, QMIN(J1), FF(J1)

            DO J2 = 1, NATOMS/2

               J3    = 3*J2
               J5    = 3*NATOMS/2 + J3
               P(:)  = QMINP(J1,J5-2:J5)

               CALL RMDRVT(P, RMI, DRMI, DRMI, DRMI, GTEST)

               DO J4 = 1, NRBSITES

                  RBCOORDS(3*J4-2:3*J4) = QMINP(J1,J3-2:J3) + MATMUL(RMI(:,:),SITE(J4,:))

               ENDDO

               DO J4 = 1, NRBSITES !- 1

                  J3 = J4 + 1
                  IF (J4 == NRBSITES) J3 = 1
!                  IF (J4 == NRBSITES - 1) J3 = 1
                  P(:) = RBCOORDS(3*J3-2:3*J3) - RBCOORDS(3*J4-2:3*J4)
                  WRITE(28,'(A4,3F20.10,2X,A12,2X,3F20.10)')                 
     &            'O', RBCOORDS(3*J4-2), RBCOORDS(3*J4-1), RBCOORDS(3*J4), 'atom_vector', P(1), P(2), P(3)

               ENDDO

            ENDDO

         ENDDO

         CLOSE(UNIT=28)

         RETURN

      ELSE IF (MULTPAHAT) THEN

         CALL VIEWMULTPAHA()
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

               CALL RMDRVT(P, RMI, DRMI, DRMI, DRMI, GTEST)

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

      ELSE IF (NTIPT) THEN

         CALL VIEWNEWTIP()
         RETURN

!|gd351>

      ELSE IF (PATCHY) THEN

         CALL VIEWPATCHY()
         RETURN

!<gd351|

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

               CALL RMDRVT(P, RMI, DRMI, DRMI, DRMI, GTEST)

               PHI = DATAN2(RMI(2,3),RMI(1,3))
               IF (PHI <= 0.D0) PHI = PHI + 2.D0*PI

               THT = DACOS(RMI(3,3))

               CHI = - DATAN2(RMI(3,2),RMI(3,1))
               IF (CHI <= 0.D0) CHI = CHI + 2.D0*PI

               PHI = PHI*180.D0/PI
               THT = THT*180.D0/PI
               CHI = CHI*180.D0/PI

               WRITE(26,'(a5,2x,3f20.10,2x,a8,6f20.10)')                 
     &              'O', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3),         
     &              'ellipse', 2.D0*PYA1(1), 2.D0*PYA1(2), 2.D0*PYA1(3), PHI, THT, CHI 

            ENDDO

         ENDDO

         CLOSE(UNIT=26)

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

               CALL RMDRVT(P, RMI, DRMI, DRMI, DRMI, GTEST)

               RBCOORDS(1:3) = MATMUL(RMI(:,:),DU(:))
               WRITE(26,'(A4,3F20.10,2X,A12,2X,3F20.10)') 'O', QMINP(J1,J3-2), QMINP(J1,J3-1), QMINP(J1,J3), 
     &         'atom_vector', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)


            ENDDO

         ENDDO

         CLOSE (UNIT=26)

         RETURN

      ELSE IF (TDHDT) THEN

         CALL VIEWTDHD()
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

               CALL RMDRVT(P, RMI, DRMI, DRMI, DRMI, GTEST)

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
               CALL RBIO(QMINP(J1,3*(J2-1)+1),QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3),
     1                       QMINP(J1,3*(NATOMS/2+J2-1)+1),
     2                       QMINP(J1,3*(NATOMS/2+J2-1)+2),
     3                       QMINP(J1,3*(NATOMS/2+J2-1)+3),
     4                       RBCOORDS,NRBSITES,SITE)
               DO J3=1,NRBSITES
                  WRITE(26,'(A4,3F20.10)') 'LA ',RBCOORDS(3*(J3-1)+1),RBCOORDS(3*(J3-1)+2),RBCOORDS(3*(J3-1)+3)
               ENDDO
            ENDDO
         ENDDO
         CLOSE(26)

      ENDIF

      RETURN
      END

      SUBROUTINE amberdump(J1,QMINP)
      USE commons
      USE modamber
      IMPLICIT NONE


      CHARACTER(LEN=25) coordfile
      CHARACTER(LEN=2) FNAME
      INTEGER J1
      DOUBLE PRECISION QMINP(NSAVE,3*NATOMS)

      IF (J1.LT.10) THEN
         WRITE (FNAME,'(I1)') J1
      ELSE
         WRITE (FNAME,'(I2)') J1
      ENDIF

      DO a=1,atoms
        x(a)=QMINP(J1,3*a-2)
        y(a)=QMINP(J1,3*a-1)
        z(a)=QMINP(J1,3*a)
      END DO

      coordfile='acoords.dump.'//FNAME

      OPEN (UNIT=4,IOSTAT=ios,FILE=coordfile,STATUS='UNKNOWN')

      DO a=1,atoms
        WRITE (UNIT=4,FMT='(A1,2X,A2,2X,I3,2X,I3,2X,F7.3,3X,F7.3,3X,F7.3)') label(a),typech(a),
     1        a,bondedto(a),x(a),y(a),z(a)
      ENDDO

      WRITE (UNIT=4,FMT='(A3)') 'end'
      WRITE (UNIT=4,FMT='(A)') ' '
      WRITE (UNIT=4,FMT='(A4,7X,I2)') 'loop',rings

      DO a=1,rings
        WRITE (UNIT=4,FMT='(I3,4X,I3)') loopatom(2*a-1),loopatom(2*a)
      END DO

      WRITE (UNIT=4,FMT='(A)') ' '
      WRITE (UNIT=4,FMT='(A7)') 'charges'

      DO a=1,atoms
        q(a)=q(a)/18.2223
        WRITE (UNIT=4,FMT='(I3,2X,F7.4)') a,q(a)
      END DO

      WRITE (UNIT=4,FMT='(A3)') 'end'

      RETURN

      END
C
C  SUBROUTINE to convert capsid CofM and DV coordinates to penatgons.
C
      SUBROUTINE CAPSIDIO(X1, Y1, Z1, L1, M1, N1,COORDS,RAD,HEIGHT)
      IMPLICIT NONE
      DOUBLE PRECISION X1, Y1, Z1, COORDS(*), HEIGHT, C2A1,
     2                 M1, L1, N1, ALPHA1, RAD, CA1, S1, C3A1,
     3                 NUM1, NUM2, NUM3, NUM4, NUM5, L12, M12, N12

      NUM1=-(1.0D0+SQRT(5.0D0))/4.0D0
      NUM2=SQRT((5.0D0-SQRT(5.0D0))/2.0D0)/2.0D0
      NUM3=SQRT((5.0D0+SQRT(5.0D0))/2.0D0)/2.0D0
      NUM4=(SQRT(5.0D0)-1.0D0)/4.0D0
      NUM5=-(1.0D0+SQRT(5.0D0))/4.0D0

      L12=L1**2
      M12=M1**2
      N12=N1**2
      ALPHA1=SQRT(L12+M12+N12)
      CA1=COS(ALPHA1)
      C2A1=RAD*CA1
      IF (ALPHA1.LT.0.0001D0) THEN
C        C3A1=RAD*(-ALPHA1/2+ALPHA1**3/24)
         C3A1=RAD*(-0.5D0+ALPHA1**2/24.0D0)
         S1=RAD*(1.0D0-ALPHA1**2/6)
      ELSE
         C3A1=RAD*(CA1-1.0D0)/ALPHA1**2
         S1=RAD*SIN(ALPHA1)/ALPHA1
      ENDIF

      COORDS(1) =     c2a1 - c3a1*l12 + x1
      COORDS(2) =     -(c3a1*l1*m1) - n1*s1 + y1
      COORDS(3) =     -(c3a1*l1*n1) + m1*s1 + z1
      COORDS(4) =     c2a1*num4 - c3a1*l1*(m1*num3 + l1*num4) + n1*num3*s1 + x1
      COORDS(5) =     c2a1*num3 - c3a1*m1*(m1*num3 + l1*num4) - n1*num4*s1 + y1
      COORDS(6) =     -(c3a1*n1*(m1*num3 + l1*num4)) - l1*num3*s1 + m1*num4*s1 + z1
      COORDS(7) =     c2a1*num1 - c3a1*l1*(l1*num1 + m1*num2) + n1*num2*s1 + x1
      COORDS(8) = c2a1*num2 - c3a1*m1*(l1*num1 + m1*num2) - n1*num5*s1 + y1
      COORDS(9) = -(c3a1*n1*(l1*num1 + m1*num2)) + m1*num1*s1 - l1*num2*s1 + z1
      COORDS(10) = c2a1*num1 + c3a1*l1*(-(l1*num1) + m1*num2) - n1*num2*s1 + x1
      COORDS(11) = -(c2a1*num2) + c3a1*m1*(-(l1*num1) + m1*num2) - n1*num5*s1 + y1
      COORDS(12) = -(c3a1*l1*n1*num1) + c3a1*m1*n1*num2 + m1*num1*s1 + l1*num2*s1 + z1
      COORDS(13) = c2a1*num4 + c3a1*l1*(m1*num3 - l1*num4) - n1*num3*s1 + x1
      COORDS(14) = -(c2a1*num3) + c3a1*m1*(m1*num3 - l1*num4) - n1*num4*s1 + y1
      COORDS(15) = c3a1*n1*(m1*num3 - l1*num4) + l1*num3*s1 + m1*num4*s1 + z1
C     COORDS(16)= (-(c3a1*l1*n1) - m1*s1 + 2*x1)/2.
C     COORDS(17)= -(c3a1*m1*n1)/2. + (l1*s1)/2. + y1
C     COORDS(18)= (c2a1 - c3a1*n12 + 2*z1)/2.
      COORDS(16)= -(c3a1*height*l1*n1) - height*m1*s1 + x1
      COORDS(17)= -(c3a1*height*m1*n1) + height*l1*s1 + y1
      COORDS(18)= c2a1*height - c3a1*height*n12 + z1

      RETURN
      END
C
C  Subroutine to convert rigid body CofM and DV coordinates to molecular sites.
C
      SUBROUTINE RBIO(X1, Y1, Z1, L1, M1, N1, COORDS, NRBSITES, SITE)
      IMPLICIT NONE
      INTEGER NRBSITES
      DOUBLE PRECISION X1, Y1, Z1, COORDS(*), C2A1, SITE(NRBSITES,3),
     2                 M1, L1, N1, ALPHA1, CA1, S1, C3A1, L12, M12, N12
      INTEGER J1

      L12=L1**2
      M12=M1**2
      N12=N1**2
      ALPHA1=SQRT(L12+M12+N12)
      CA1=COS(ALPHA1)
      C2A1=CA1
      IF (ALPHA1.LT.0.0001D0) THEN
C        C3A1=(-ALPHA1/2+ALPHA1**3/24)
         C3A1=(-0.5D0+ALPHA1**2/24.0D0)
         S1=(1.0D0-ALPHA1**2/6)
      ELSE
         C3A1=(CA1-1.0D0)/ALPHA1**2
         S1=SIN(ALPHA1)/ALPHA1
      ENDIF
   
      DO J1=1,NRBSITES
         COORDS(3*(J1-1)+1)=c2a1*SITE(J1,1) + s1*(n1*SITE(J1,2) - m1*SITE(J1,3)) - 
     1                      c3a1*l1*(l1*SITE(J1,1) + m1*SITE(J1,2) + n1*SITE(J1,3)) + X1
         COORDS(3*(J1-1)+2)=c2a1*SITE(J1,2) + s1*(-(n1*SITE(J1,1)) + l1*SITE(J1,3)) 
     1                    - c3a1*m1*(l1*SITE(J1,1) + m1*SITE(J1,2) + n1*SITE(J1,3)) + Y1
         COORDS(3*(J1-1)+3)=s1*(m1*SITE(J1,1) - l1*SITE(J1,2)) + c2a1*SITE(J1,3) 
     1                    - c3a1*n1*(l1*SITE(J1,1) + m1*SITE(J1,2) + n1*SITE(J1,3)) + Z1
      ENDDO 

      RETURN
      END
C
C  SUBROUTINE to convert TIP oxygen and DV coordinates to Cartesians.
C
      SUBROUTINE TIPIO(X1, Y1, Z1, L1, M1, N1, COORDS)
      IMPLICIT NONE
      DOUBLE PRECISION X1, Y1, Z1, COORDS(*), C2A1, M1, L1, N1, ALPHA1, CA1, S1, C3A1, L12, M12, N12

      L12=L1**2
      M12=M1**2
      N12=N1**2
      ALPHA1=SQRT(L12+M12+N12)
      CA1=COS(ALPHA1)
      C2A1=CA1
      IF (ALPHA1.LT.0.0001D0) THEN
C        C3A1=(-ALPHA1/2+ALPHA1**3/24)
         C3A1=(-0.5D0+ALPHA1**2/24.0D0)
         S1=(1.0D0-ALPHA1**2/6)
      ELSE
         C3A1=(CA1-1.0D0)/ALPHA1**2
         S1=SIN(ALPHA1)/ALPHA1
      ENDIF

      COORDS(1) = X1
      COORDS(2) = Y1
      COORDS(3) = Z1    
      COORDS(4) = 0.756950327*c2a1 - c3a1*l1*(0.756950327*l1 - 0.585882276*n1) + 0.585882276*m1*s1 + X1
      COORDS(5) = -(c3a1*m1*(0.756950327*l1 - 0.585882276*n1)) + (-0.585882276*l1 - 0.756950327*n1)*s1 + Y1
      COORDS(6) = -0.585882276*c2a1 - c3a1*(0.756950327*l1 - 0.585882276*n1)*n1 + 0.756950327*m1*s1 + Z1
      COORDS(7) = -0.756950327*c2a1 + c3a1*l1*(0.756950327*l1 + 0.585882276*n1) + 0.585882276*m1*s1 + X1
      COORDS(8) = c3a1*m1*(0.756950327*l1 + 0.585882276*n1) + (-0.585882276*l1 + 0.756950327*n1)*s1 + Y1
      COORDS(9) = -0.585882276*c2a1 + c3a1*(0.756950327*l1 + 0.585882276*n1)*n1 - 0.756950327*m1*s1 + Z1

      RETURN
      END

