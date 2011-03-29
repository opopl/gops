!OP226>=================================== 
!OP226> GPL LICENSE INFO {{{ 
C   GMIN: A PROGRAM FOR FINDING GLOBAL MINIMA
C   COPYRIGHT (C) 1999-2006 DAVID J. WALES
C   THIS FILE IS PART OF GMIN.
C
C   GMIN IS FREE SOFTWARE; YOU CAN REDISTRIBUTE IT AND/OR MODIFY
C   IT UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE AS PUBLISHED BY
C   THE FREE SOFTWARE FOUNDATION; EITHER VERSION 2 OF THE LICENSE, OR
C   (AT YOUR OPTION) ANY LATER VERSION.
C
C   GMIN IS DISTRIBUTED IN THE HOPE THAT IT WILL BE USEFUL,
C   BUT WITHOUT ANY WARRANTY; WITHOUT EVEN THE IMPLIED WARRANTY OF
C   MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  SEE THE
C   GNU GENERAL PUBLIC LICENSE FOR MORE DETAILS.
C
C   YOU SHOULD HAVE RECEIVED A COPY OF THE GNU GENERAL PUBLIC LICENSE
C   ALONG WITH THIS PROGRAM; IF NOT, WRITE TO THE FREE SOFTWARE
C   FOUNDATION, INC., 59 TEMPLE PLACE, SUITE 330, BOSTON, MA  02111-1307  USA
C
!OP226>}}} 
!OP226>=================================== 
      SUBROUTINE FINALIO
!OP226>=================================== 
!> \BRIEF PRODUCE FINAL QUENCHES  
!OP226>=================================== 
!OP226> DECLARATIONS {{{ 
      USE COMMONS
      USE MODAMBER
      USE MODAMBER9, ONLY : COORDS1,LCRD,IH,M04,NATOM,AMBFINALIO_NODE
      USE PYMODULE, ONLY : SITECOORDS,ELLST1,ELLMAT
      USE QMODULE
      USE MODCHARMM
      USE AMHGLOBALS, ONLY:NMRES,IRES

      IMPLICIT NONE

C   MCP
      INTEGER III, I3,  GLY_COUNT, ID, NUMCRD, NUMPRO, NCPHST
      INTEGER J1, J2, J3, J4, J5, NTYPEA, MYUNIT2, I1, NDUMMY, MYUNIT3, NC, NRBS1, NRBS2
      DOUBLE PRECISION EPSAB, EPSBB, SIGAB, SIGBB, RBCOORDS(NRBSITES*3), DCOORDS(3*NATOMS)
      DOUBLE PRECISION P3(3,3), P(3), DU(3), RMI(3,3), DRMI(3,3), PI, PHI, THT, CHI
      DOUBLE PRECISION, ALLOCATABLE :: XCOORDS(:), YCOORDS(:)
      COMMON /BIN/ EPSAB, EPSBB, SIGAB, SIGBB, NTYPEA
      CHARACTER(LEN=13) J1CHAR,J1CHAR2                  !FOR GAY-BERNE OUTPUT FILES
      DOUBLE PRECISION EULERPHI,EULERPSI,EULERTHETA,EULERTHETADEG,EULERPHIDEG,EULERPSIDEG  ! EULER ANGLES FOR ELLIPSOIDS OF REVOLUTION

      DOUBLE PRECISION EPS2, RAD, HEIGHT,SUMX,SUMY,SUMZ
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
!OP226> END DECLARATIONS }}} 

      PI = 4.D0*DATAN(1.D0)

      NUMPRO = 1
      NUMCRD = 3

      ALLOCATE(DBNAME(NSAVE))

     
      DO J1=1,NSAVE
         WRITE(DBNUM,*) J1
         DBNAME(J1)='DBASE.'//TRIM(ADJUSTL(DBNUM))
      ENDDO
      DCOORDS(1:3*NATOMS)=0.0D0
!      IF (AMH) THEN
!         CALL WALESAMH_FINALIO
!         STOP
!      ENDIF

       IF (AMHT) THEN
         OPEN(UNIT=26,FILE='MOVIE_GMIN',STATUS='UNKNOWN')
         WRITE(26,334)NMRES,NUMCRD,NUMPRO,NSAVE
334      FORMAT(4(I8,1X),' NMRES NMCRD NUMPRO NMSNAP')
       ENDIF


      IF (MPIT) THEN
         WRITE (ISTR, '(I10)') MYUNIT-22980+1 
         MYUNIT2=(MYUNIT-22980+1)+100
         MYFILENAME2="LOWEST."//TRIM(ADJUSTL(ISTR))
         OPEN(MYUNIT2,FILE=TRIM(ADJUSTL(MYFILENAME2)), STATUS="UNKNOWN", FORM="FORMATTED")
         IF (CHRMMT) THEN
            DO J1=1,NSAVE
               WRITE(DBNUM,*) J1
               DBNAME(J1)='DBASE.'//TRIM(ADJUSTL(ISTR))//'.'//TRIM(ADJUSTL(DBNUM))             
            ENDDO
         ENDIF
         IF (CSMT.AND.(.NOT.SYMMETRIZECSM)) THEN
            MYUNIT3=(MYUNIT+22980+2)+111
            MYFILENAME3="CSMAV."//TRIM(ADJUSTL(ISTR))
         ENDIF
      ELSE
         MYUNIT2=25 
         OPEN(MYUNIT2,FILE='LOWEST',STATUS='UNKNOWN')
         IF (CSMT.AND.(.NOT.SYMMETRIZECSM)) THEN
            MYUNIT3=26 
            OPEN(MYUNIT3,FILE='CSMAV.XYZ',STATUS='UNKNOWN')
         ENDIF
      ENDIF

      DO J1=1,NSAVE
         IF (AMHT) THEN
            COUNTT=J1
            CALL NUM_TO_CHAR(COUNTT,COUNTTT)
            OPEN(UNIT=27,FILE='MOVIE_GMIN.'//COUNTTT//'.PDB',STATUS='UNKNOWN')
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
!        IF (CSMT.AND.DEBUG) WRITE(MYUNIT,'(A,I6,2G20.10)') 'FINALIO> J1,QMIN,QMINAV=',J1,QMIN(J1),QMINAV(J1)
         WRITE(MYUNIT2,10) J1, QMIN(J1), FF(J1)
10       FORMAT('ENERGY OF MINIMUM ',I6,'=',F20.10,' FIRST FOUND AT STEP ',I8)
         IF (MSORIGT.OR.FRAUSIT) THEN
            WRITE(MYUNIT2,20) (QMINP(J1,J2),J2=1,3*(NATOMS-NS))
20          FORMAT('SI',3F20.10)
         ELSE IF (MSTRANST) THEN
            WRITE(MYUNIT2,20) (QMINP(J1,J2),J2=1,3*(NATOMS-NS))
         ELSE IF (RGCL2) THEN
            WRITE(MYUNIT2,'(A,F20.10)') 'CL 0.0 0.0 ', 0.995D0
            WRITE(MYUNIT2,'(A,F20.10)') 'CL 0.0 0.0 ',-0.995D0
            WRITE(MYUNIT2,60) (QMINP(J1,J2),J2=1,3*(NATOMS-NS))
60          FORMAT('AR ',3F20.10)
         ELSE IF (AMHT) THEN
C
C   OUTPUT CORDS FOR LOWEST IN X,Y,Z FORMAT
C   OUTPUT CORDS FOR MOVIE_GMIN IN MOVIESEG FORMAT
C
            WRITE(26,683)NUMCRD,J1,NUMCRD,REAL(NUMCRD),NUMCRD
683         FORMAT(3(I6,1X),F8.4,1X,I5,' STUCT SNAP T T TID')

            GLY_COUNT = 0

            DO 1964 III = 1,NMRES
               IF (IRES(III).EQ.8) THEN
!!                PPPCORD(RESIDUE, XYZ, NUMPRO, ATOM TYPES
                  PPPCORD(III, 1, 1, 1) = REAL(QMINP(J1,9*(III-1)+1- GLY_COUNT*3)) !  CA X
                  PPPCORD(III, 2, 1, 1) = REAL(QMINP(J1,9*(III-1)+2- GLY_COUNT*3)) !  CA Y
                  PPPCORD(III, 3, 1, 1) = REAL(QMINP(J1,9*(III-1)+3- GLY_COUNT*3)) !  CA Z
!    SWAP  CA FOR CB 
                  PPPCORD(III, 1, 1, 2) = REAL(QMINP(J1,9*(III-1)+1- GLY_COUNT*3)) !  CB X
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
                  PPPCORD(III, 1, 1, 1) = REAL(QMINP(J1,9*(III-1)+1- GLY_COUNT*3)) !  CA X
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
!632           FORMAT('CA: ',3(F25.15,1X),'CB: ',3(F25.15,1X),'OX: ', 3(F25.15,1X))
526         CONTINUE

            DO III = 1+1, NMRES
               PPPCORD(III,1,1,4) =
     &          0.4831806D0*PPPCORD(III-1,1,1,1) + 0.7032820D0*PPPCORD(III,1,1,1) - 0.1864626D0*PPPCORD(III-1,1,1,3)
               PPPCORD(III,2,1,4) =
     &          0.4831806D0*PPPCORD(III-1,2,1,1) + 0.7032820D0*PPPCORD(III,2,1,1) - 0.1864626D0*PPPCORD(III-1,2,1,3)
               PPPCORD(III,3,1,4) =
     &          0.4831806D0*PPPCORD(III-1,3,1,1) + 0.7032820D0*PPPCORD(III,3,1,1) - 0.1864626D0*PPPCORD(III-1,3,1,3)
            ENDDO

            DO III = 1, NMRES-1
               PPPCORD(III,1,1,5) =
     &          0.4436538D0*PPPCORD(III,1,1,1)+0.2352006D0*PPPCORD(III+1,1,1,1)+0.3211455D0*PPPCORD(III,1,1,3)
               PPPCORD(III,2,1,5) =
     &          0.4436538D0*PPPCORD(III,2,1,1)+0.2352006D0*PPPCORD(III+1,2,1,1)+0.3211455D0*PPPCORD(III,2,1,3)
               PPPCORD(III,3,1,5) =
     &          0.4436538D0*PPPCORD(III,3,1,1)+0.2352006D0*PPPCORD(III+1,3,1,1)+0.3211455D0*PPPCORD(III,3,1,3)
            ENDDO

            DO III = 1, NMRES
               RES_TYPE = AMINOA(IRES(III))

               IF (III .NE. 1) THEN
                  ATOM_TYPE='N '
                  ID = 4
                  WRITE(27,61)III,ATOM_TYPE,RES_TYPE,III,PPPCORD(III,1,1,ID),PPPCORD(III,2,1,ID),PPPCORD(III,3,1,ID),III
61      FORMAT('ATOM',4X,I3,2X,A2,2X,A3,3X,I3,4X,F8.3,F8.3,F8.3,2X,'1.00',2X,'0.00',6X,'TPDB',1X,I3)

               ENDIF

               ATOM_TYPE='CA'
               ID = 1
               WRITE(27,61)III,ATOM_TYPE,RES_TYPE,III,PPPCORD(III,1,1,ID),PPPCORD(III,2,1,ID),PPPCORD(III,3,1,ID),III

               IF (RES_TYPE .NE. 'GLY') THEN
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
               IF (ZSYM(J2).EQ.'PL') WRITE(MYUNIT2,'(A,3F20.10)') 'NA  ',(QMINP(J1,3*(J2-1)+J3),J3=1,3)
               IF (ZSYM(J2).EQ.'MI') WRITE(MYUNIT2,'(A,3F20.10)') 'CL  ',(QMINP(J1,3*(J2-1)+J3),J3=1,3)
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
               WRITE(MYUNIT2,'(A,3F20.10)') TYPECH(J2)(1:1),(QMINP(J1,3*(J2-1)+J3),J3=1,3)
            ENDDO
         ELSE IF (AMBERT) THEN
          ! SF344> WRITE OUT COORDINATES
            COORDS1(1:3*NATOMS) = QMINP(J1,1:3*NATOMS)
            IF (DUMPSTRUCTURES) THEN
              CALL AMBERFINALIO(J1,25,AMBFINALIO_NODE,'0',0,COORDS1(1:3*NATOMS))         
              WRITE(J1CHAR2,'(I3)') J1
              WRITE(J1CHAR,'(A,A)') 'COORDS.',TRIM(ADJUSTL(J1CHAR2))
              OPEN(UNIT=226,FILE=TRIM(ADJUSTL(J1CHAR)),STATUS='UNKNOWN')

              DO J2=1,NATOMS
                 WRITE(226,'(3F28.20)') QMINP(J1,3*(J2-1)+1),QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3)
              ENDDO
              CLOSE(226)
               WRITE(J1CHAR2,'(I3)') J1
               WRITE(J1CHAR,'(A,A)') 'COORDS.',TRIM(ADJUSTL(J1CHAR2))
               OPEN(UNIT=226,FILE=TRIM(ADJUSTL(J1CHAR)),STATUS='UNKNOWN')
              
               DO J2=1,NATOMS
                  WRITE(226,'(3F28.20)') QMINP(J1,3*(J2-1)+1),QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3)
               ENDDO
               CLOSE(226)

            ELSE
               DO I1=1,NATOMS
                  WRITE(MYUNIT2,'(A2,3F20.10)') IH(M04+I1-1),COORDS1(3*I1-2),COORDS1(3*I1-1),COORDS1(3*I1)
               ENDDO
            ENDIF

         ELSE IF (ELLIPSOIDT.OR.LJCAPSIDT.OR.GBT.OR.GBDT) THEN
C
C DETERMINE THE CENTRE OF COORDINATES, THEN CENTRE THEM
C
              SUMX=0.0D0
              SUMY=0.0D0
              SUMZ=0.0D0

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

              DO J2=1,NATOMS/2
                 IF (PARAMONOVPBCX) THEN
                        ! ENSURE Y COMPONENT OF PARTICLE 1 VECTOR IS WITHIN BOXLY/2 OF ZERO. 
                        ! IF IT ISN'T THEN SUBTRACT INTEGER NUMBER OF BOXLY'S SUCH THAT IT IS.
                    QMINP(J1,3*J2-2)=QMINP(J1,3*J2-2)-BOXLX*NINT(QMINP(J1,3*J2-2)/BOXLX)
                 ENDIF
                 IF (PARAMONOVPBCY) THEN
                        ! ENSURE Y COMPONENT OF PARTICLE 1 VECTOR IS WITHIN BOXLY/2 OF ZERO. 
                        ! IF IT ISN'T THEN SUBTRACT INTEGER NUMBER OF BOXLY'S SUCH THAT IT IS.
                    QMINP(J1,3*J2-1)=QMINP(J1,3*J2-1)-BOXLY*NINT(QMINP(J1,3*J2-1)/BOXLY)
                 ENDIF
                 IF (PARAMONOVPBCZ) THEN
                        ! ENSURE Y COMPONENT OF PARTICLE 1 VECTOR IS WITHIN BOXLY/2 OF ZERO. 
                        ! IF IT ISN'T THEN SUBTRACT INTEGER NUMBER OF BOXLY'S SUCH THAT IT IS.
                    QMINP(J1,3*J2  )=QMINP(J1,3*J2  )-BOXLZ*NINT(QMINP(J1,3*J2  )/BOXLZ)
                 ENDIF
              ENDDO

              DO J2=1,NATOMS/2
                 WRITE(MYUNIT2,'(A5,2X,3F20.10,2X,A11,3F20.10)') 'H',QMINP(J1,3*(J2-1)+1), 
     &             QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3),
     &             'ATOM_VECTOR',QMINP(J1,3*NATOMS/2+3*(J2-1)+1),
     &             QMINP(J1,3*NATOMS/2+3*(J2-1)+2),QMINP(J1,3*NATOMS/2+3*(J2-1)+3)
             ENDDO

         ELSE IF (CHRMMT) THEN
            DO J2=1,NATOMS
               WRITE(MYUNIT2,'(A,1X,3F20.10)') ZSYM(J2)(1:1),(QMINP(J1,3*(J2-1)+J3),J3=1,3)
            ENDDO
!       CSW34> THIS DO LOOP APPEARED TO BE BE MISSING ON 30/9/08 WHICH WOULD EASILY
!       EXPLAIN THE OUTPUT PROBLEMS! 
            DO J2=1,NATOMS
               DCOORDS(3*(J2-1)+1)=QMINP(J1,3*(J2-1)+1)
               DCOORDS(3*(J2-1)+2)=QMINP(J1,3*(J2-1)+2)
               DCOORDS(3*(J2-1)+3)=QMINP(J1,3*(J2-1)+3)
            ENDDO
            CALL CHARMMDUMP(DCOORDS,DBNAME(J1))

!    DC430 >
!    |GD351> ADDED PATCHY

         ELSE IF (DBPT .OR. DBPTDT .OR. DMBLMT .OR. LINRODT .OR. LWOTPT .OR. MSTBINT .OR. MSSTOCKT .OR. NCAPT .OR. NPAHT  
     &           .OR. NTIPT .OR. STOCKAAT .OR. PAHAT .OR. PAHW99T .OR. TDHDT .OR. WATERDCT .OR. WATERKZT .OR. PATCHY) THEN
            DO J2 = 1, NATOMS/2
               WRITE(MYUNIT2,'(3F20.10)') (QMINP(J1,3*(J2-1)+J3),J3=1,3)
            ENDDO
            DO J2 = 1, NATOMS/2
               WRITE(MYUNIT2,'(3F20.10)') (QMINP(J1,3*NATOMS/2+3*(J2-1)+J3),J3=1,3)
            ENDDO

         ELSE IF (GBT.OR.GBDT.OR.GBDPT.OR.PYGT.OR.PYGDPT.OR.MSGBT.OR.MSPYGT) THEN
            DO J2 = 1, NATOMS/2
               WRITE(MYUNIT2,'(3F20.10)') (QMINP(J1,3*(J2-1)+J3),J3=1,3)    
            ENDDO
            DO J2 = 1, NATOMS/2
               WRITE(MYUNIT2,'(3F20.10)') (QMINP(J1,3*NATOMS/2+3*(J2-1)+J3),J3=1,3)
            ENDDO

         ELSE IF (GEMT) THEN
            DO J2 = 1, NATOMS
               WRITE(MYUNIT2,'(3F20.10)') (QMINP(J1,3*(J2-1)+J3),J3=1,3)
            ENDDO

         ELSE IF (BLNT.AND.(.NOT.P46).AND.(.NOT.G46)) THEN
C
C THIS WRITES 'LOWEST' IN XYZ (XMAKEMOL) FORMAT
C
            DO J2=1,NATOMS
               WRITE(MYUNIT2,'(2A1,1X,3F20.10)') BEADLETTER(J2),'L',(QMINP(J1,3*(J2-1)+J3),J3=1,3)
            ENDDO
         ELSE
            IF (CSMT.AND.(.NOT.SYMMETRIZECSM)) THEN
               WRITE(MYUNIT3,'(I6)') NATOMS
               WRITE(MYUNIT3,'(A,I6,2(A,G20.10))') 'AVERAGED STRUCTURE FOR FINAL SOLUTION ',J1,
     &                                          ' CSM=',QMINAV(J1),' CSM FOR REFERENCE STRUCTURE=',QMIN(J1)
               WRITE(MYUNIT3,30) (QMINPCSMAV(J1,J2),J2=1,3*(NATOMS-NS))
            ENDIF
            WRITE(MYUNIT2,30) (QMINP(J1,J2),J2=1,3*(NATOMS-NS))
30          FORMAT('LA ',3F20.10)
         ENDIF

!|GD351>
         IF (ASAOOS) THEN
            OPEN(31,FILE='PARTICLES.XYZ')
            WRITE(31,*) NATOMS
            WRITE(31,*) ' '
            DO J2=1,NATOMS
                  WRITE(31,'(A,3F20.10)') 'H  ',(QMINP(J1,3*(J2-1)+J3),J3=1,3)
            ENDDO
            CLOSE(31)
         END IF
!<GD351|

         IF ((NS.GT.0).AND.(.NOT.(WELCH.OR.TOSI))) THEN
            IF (MSORIGT.OR.FRAUSIT) THEN
               WRITE(MYUNIT2,40) (QMINP(J1,J2),J2=3*(NATOMS-NS)+1,3*NATOMS)
40             FORMAT('SI',3F20.10)
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
C  END OF LOOP OVER DUMP TO FILE LOWEST OR EQUIVALENT.
C
      CLOSE(MYUNIT2)
      IF (CSMT.AND.(.NOT.SYMMETRIZECSM)) CLOSE(MYUNIT3)
C
C     CSW34> NEW LOOP FOR DUMPING INTERACTION ENERGY FILES IF A9INTE IS SPECIFIED
C     ADDED THE MISSING IF BLOCK TO TEST FOR A9INTE 9/12/09 DJW
C
      IF (A9INTET) THEN
         IF (MPIT) THEN
            WRITE (ISTR, '(I10)') MYUNIT-22980+1
            MYUNIT2=(MYUNIT-22980+1)+100
            MYFILENAME2="INTELOWEST."//TRIM(ADJUSTL(ISTR))
            OPEN(MYUNIT2,FILE=TRIM(ADJUSTL(MYFILENAME2)), STATUS="UNKNOWN", FORM="FORMATTED")
         ELSE
            MYUNIT2=25
            OPEN(MYUNIT2,FILE='INTELOWEST',STATUS='UNKNOWN')
         ENDIF
      ENDIF
C
C     CSW34> LOOP STRUCTURE COPIED FROM THE ELSEIF(AMBERT) BLOCK ABOVE
C
      IF (A9INTET.AND.AMBERT) THEN
         DO J1=1,NSAVEINTE
            WRITE(MYUNIT2,*) NATOMS
!     CSW34> WRITE HEADER TO INTELOWEST FOR CURRENT MINIMUM
            WRITE(MYUNIT2,10) J1, INTEQMIN(J1), INTEFF(J1)
!     SF344> WRITE OUT COORDINATES
            COORDS1(1:3*NATOMS) = INTEQMINP(J1,1:3*NATOMS)
            IF (DUMPSTRUCTURES) THEN
              CALL INTEFINALIO(J1,MYUNIT2,AMBFINALIO_NODE,'0',0,COORDS1(1:3*NATOMS))         
              WRITE(J1CHAR2,'(I3)') J1
              WRITE(J1CHAR,'(A,A)') 'INTECOORDS.',TRIM(ADJUSTL(J1CHAR2))
              OPEN(UNIT=226,FILE=TRIM(ADJUSTL(J1CHAR)),STATUS='UNKNOWN')

              DO J2=1,NATOMS
                 WRITE(226,'(3F28.20)') INTEQMINP(J1,3*(J2-1)+1),INTEQMINP(J1,3*(J2-1)+2),INTEQMINP(J1,3*(J2-1)+3)
              ENDDO
              CLOSE(226)
               WRITE(J1CHAR2,'(I3)') J1
               WRITE(J1CHAR,'(A,A)') 'INTECOORDS.',TRIM(ADJUSTL(J1CHAR2))
               OPEN(UNIT=226,FILE=TRIM(ADJUSTL(J1CHAR)),STATUS='UNKNOWN')
              
               DO J2=1,NATOMS
                  WRITE(226,'(3F28.20)') INTEQMINP(J1,3*(J2-1)+1),INTEQMINP(J1,3*(J2-1)+2),INTEQMINP(J1,3*(J2-1)+3)
               ENDDO
               CLOSE(226)

            ELSE
               DO I1=1,NATOMS
                  WRITE(MYUNIT2,'(A2,3F20.10)') IH(M04+I1-1),COORDS1(3*I1-2),COORDS1(3*I1-1),COORDS1(3*I1)
               ENDDO
            ENDIF
         ENDDO
      ENDIF
C
C  END OF LOOP OVER DUMP TO FILE INTELOWEST 
C
      CLOSE(MYUNIT2)

!     CSW34> EDITS TO THE RMS KEYWORD          
      IF (CHRMMT.AND.RMST) THEN
!        IF (PROGRESS) THEN
!           DCOORDS(1:3*NATOMS)=RMSCOOR(1,1:3*NATOMS)
!           IF(RMSBEST(1,2)<0.D0) CALL CHARMMDUMP(DCOORDS,'CLOSESTRMS')
!           WRITE(MYUNIT,'(A9,F8.5)') 'RMSDMIN= ',RMSBEST(1,1)
!        ELSE      
! REMEMBER TO RE-INDENT THE BELOW IF UNCOMMENTING ABOVE!
         OPEN(UNIT=MYUNIT2,FILE='RMSBEST',STATUS='UNKNOWN')
         DO J2=1,RMSSAVE
            WRITE(MYUNIT2,'(I6,F6.3,F15.5)')J2,RMSBEST(J2,1),RMSBEST(J2,2)
            WRITE(CRMS,'(I6)') J2
            DCOORDS(1:3*NATOMS)=RMSCOOR(J2,1:3*NATOMS)
            IF(RMSBEST(J2,2)<0.D0) CALL CHARMMDUMP(DCOORDS,'RMS.'//TRIM(ADJUSTL(CRMS)))
         ENDDO
         CLOSE(MYUNIT2)
!        ENDIF
      ENDIF

      IF (LJCOULT) THEN
         OPEN(UNIT=26,FILE='LJCOUL.XYZ',STATUS='UNKNOWN')
         DO J1=1,NSAVE
            WRITE(26,'(I6)') NATOMS
            WRITE(26,10) J1, QMIN(J1), FF(J1)
            DO J2=1,NATOMS
C              USE "O" ATOM TYPE TO HIGHLIGHT CHARGED PARTICLES AND "N" FOR THE NEUTRAL ONES.
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
         OPEN(UNIT=26,FILE='STOCK.XYZ',STATUS='UNKNOWN')
         DO J1=1,NSAVE
            WRITE(26,'(I6)') NATOMS/2
            WRITE(26,10) J1, QMIN(J1), FF(J1)
            DO J2=1,NATOMS/2
               WRITE (26,'(A4,3F18.10,A12,3F18.10)') 'LA ',
     &            QMINP(J1,3*(J2-1)+1),QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3),
     &            'ATOM_VECTOR',
     &            SIN(QMINP(J1,3*(NATOMS/2)+3*(J2-1)+1))*COS(QMINP(J1,3*(NATOMS/2)+3*(J2-1)+2)),
     &            SIN(QMINP(J1,3*(NATOMS/2)+3*(J2-1)+1))*SIN(QMINP(J1,3*(NATOMS/2)+3*(J2-1)+2)),
     &                    COS(QMINP(J1,3*(NATOMS/2)+3*(J2-1)+1))
            ENDDO
         ENDDO
         CLOSE(26)

      ELSE IF (ELLIPSOIDT.OR.LJCAPSIDT.OR.GBT.OR.GBDT.OR.MULTISITEPYT) THEN
         DO J1=1,NSAVE
            WRITE(J1CHAR2,'(I3)') J1
!            WRITE(J1CHAR,'(A5,A,A4)') 'GBMIN',TRIM(ADJUSTL(J1CHAR2)),'.XYZ'
!            OPEN(UNIT=26,FILE=TRIM(ADJUSTL(J1CHAR)),STATUS='UNKNOWN')
!
!            DO J2=1,NATOMS/2
!               WRITE(26,'(3F18.10)') QMINP(J1,3*(J2-1)+1),QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3)
!            ENDDO
          IF (MPIT) THEN
              WRITE (ISTR, '(I10)') MYUNIT-22980+1
              MYUNIT2=(MYUNIT-22980+1)+100
              MYFILENAME2='COORDS.'//TRIM(ADJUSTL(J1CHAR2))//'.'//TRIM(ADJUSTL(ISTR))
              OPEN(MYUNIT2,FILE=TRIM(ADJUSTL(MYFILENAME2)), STATUS="UNKNOWN", FORM="FORMATTED")
          ELSE 
              MYUNIT2=226
              MYFILENAME2='COORDS.'//TRIM(ADJUSTL(J1CHAR2))
              OPEN(MYUNIT2,FILE=TRIM(ADJUSTL(MYFILENAME2)), STATUS="UNKNOWN", FORM="FORMATTED")
          END IF

            DO J2=1,NATOMS
               WRITE(MYUNIT2,'(3F28.20)') QMINP(J1,3*(J2-1)+1),QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3)
            ENDDO

            CLOSE(MYUNIT2)

         ENDDO
C
C  WRITE OUT LOWEST NSAVE STRUCTURES TO XMAKEMOL ELLIPSOID FORMAT FOR 
C  CLUSTERS OF ELLIPSOIDS OF REVOLUTION 
C
C         OPEN(UNIT=26,FILE="ELLIPSOID.XYZ",STATUS='UNKNOWN')
      IF (MPIT) THEN
         WRITE (ISTR, '(I10)') MYUNIT-22980+1
         MYUNIT2=(MYUNIT-22980+1)+300
         MYFILENAME2="ELLIPSOID."//TRIM(ADJUSTL(ISTR))//".XYZ"
         OPEN(MYUNIT2,FILE=TRIM(ADJUSTL(MYFILENAME2)), STATUS="UNKNOWN", FORM="FORMATTED")
      ELSE
         MYUNIT2=26
         OPEN(MYUNIT2,FILE='ELLIPSOID.XYZ',STATUS='UNKNOWN')
      ENDIF

        DO J1=1,NSAVE
        IF(.NOT.MULTISITEPYT) THEN
         WRITE(MYUNIT2,*) NATOMS/2
        ELSE
         WRITE(MYUNIT2,*) NATOMS*NPYSITE/2
        END IF
         WRITE(MYUNIT2,10) J1, QMIN(J1), FF(J1)

            DO J2=1,NATOMS/2
               
               IF (GAYBERNET) THEN
                  CALL ELLIPSOIDSAATOPOLAR(QMINP(J1,3*NATOMS/2+3*(J2-1)+1),QMINP(J1,3*NATOMS/2+3*(J2-1)+2),
     &                    QMINP(J1,3*NATOMS/2+3*(J2-1)+3),EULERPHI,EULERPSI,EULERTHETA,EULERPHIDEG,EULERPSIDEG,EULERTHETADEG)
!                 EULERPHIDEG = 90-EULERPHIDEG    ! EULERPHIDEG RETURNED FROM ELLIPSOIDSAATOPOLAR IS IN FACT THE ALPHA ANGLE
!                                                ! DEFINED BY ME (ANGLE OF VECTOR WITH THE XY PLANE)
                  WRITE(MYUNIT2,'(A5,2X,3F20.10,2X,A8,6F15.8,2X,A11,3F15.8)') '0',QMINP(J1,3*(J2-1)+1),
     &                  QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3),
     &                  'ELLIPSE ',1.0D0,1.0D0,GBANISOTROPYR,
     &                                                 EULERPSIDEG,EULERPHIDEG,0.0D0, ! THIS IS IN DEGREES
     &                  'ATOM_VECTOR',QMINP(J1,3*NATOMS/2+3*(J2-1)+1),
     &                  QMINP(J1,3*NATOMS/2+3*(J2-1)+2),QMINP(J1,3*NATOMS/2+3*(J2-1)+3)
               ELSE IF (PARAMONOVT) THEN
                  CALL ELLIPSOIDSAATOPOLAR(QMINP(J1,3*NATOMS/2+3*(J2-1)+1),QMINP(J1,3*NATOMS/2+3*(J2-1)+2),
     &                    QMINP(J1,3*NATOMS/2+3*(J2-1)+3),EULERPHI,EULERPSI,EULERTHETA,EULERPHIDEG,EULERPSIDEG,EULERTHETADEG)
                  WRITE(MYUNIT2,'(A5,2X,3F20.10,2X,A8,6F15.8,2X,A11,3F15.8)') 'O',QMINP(J1,3*(J2-1)+1),
     &             QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3),
     &             'ELLIPSE ',PARAMC1*2.0D0,PARAMB1*2.0D0,PARAMA1*2.0D0,EULERPHIDEG,EULERPSIDEG,EULERTHETADEG,
     &             'ATOM_VECTOR',QMINP(J1,3*NATOMS/2+3*(J2-1)+1),
     &             QMINP(J1,3*NATOMS/2+3*(J2-1)+2),QMINP(J1,3*NATOMS/2+3*(J2-1)+3)
               ELSE IF (PYGPERIODICT.OR.PYBINARYT) THEN
                  CALL AATOEULER(QMINP(J1,3*NATOMS/2+3*(J2-1)+1),QMINP(J1,3*NATOMS/2+3*(J2-1)+2),
     &                           QMINP(J1,3*NATOMS/2+3*(J2-1)+3),EULERPHIDEG,EULERPSIDEG,EULERTHETADEG)

                  WRITE(MYUNIT2,'(A5,2X,3F20.10,2X,A8,6F15.8,2X,A11,3F15.8)') 'O',QMINP(J1,3*(J2-1)+1),
     &                  QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3),
     &            'ELLIPSE ',PYA1BIN(J2,1)*2.0D0,PYA1BIN(J2,2)*2.0D0,PYA1BIN(J2,3)*2.0D0,EULERPHIDEG,EULERPSIDEG,EULERTHETADEG,
     &            'ATOM_VECTOR',QMINP(J1,3*NATOMS/2+3*(J2-1)+1),
     &            QMINP(J1,3*NATOMS/2+3*(J2-1)+2),QMINP(J1,3*NATOMS/2+3*(J2-1)+3)
               ELSE IF (MULTISITEPYT) THEN
                 CALL AATOSITES(QMINP(J1,3*(J2-1)+1:3*(J2-1)+3),QMINP(J1,3*NATOMS/2+3*(J2-1)+1:3*NATOMS/2+3*(J2-1)+3),SITECOORDS)

!                  CALL AATOEULER(QMINPS(J1,3*NPYSITE*NATOMS/2+3*(J2-1)+1),QMINP(J1,3*NPYSITE*NATOMS/2+3*(J2-1)+2),
!     &                           QMINP(J1,3*NPYSITE*NATOMS/2+3*(J2-1)+3),EULERPHIDEG,EULERPSIDEG,EULERTHETADEG)
                DO J3=1,NPYSITE
!                   CALL AATOEULER(ELLST3(J3,1),ELLST3(J3,2),ELLST3(J3,3),EULERPHIDEG,EULERTHETADEG,EULERPSIDEG)

!                  WRITE(MYUNIT2,'(A5,2X,3F20.10)') 'LA',SITECOORDS(J3,1), SITECOORDS(J3,2), SITECOORDS(J3,3)
                  WRITE(MYUNIT2,'(A5,2X,3F20.10,2X,A8,12F15.8,2X,A11,3F15.8)') 'O',
     &                                          SITECOORDS(J3,1),SITECOORDS(J3,2),SITECOORDS(J3,3),
     &            'ELLIPSE ',ELLST1(J3,1)*2.0D0,ELLST1(J3,2)*2.0D0,ELLST1(J3,3)*2.0D0,
     &                       ELLMAT(J3,1,1),ELLMAT(J3,1,2),ELLMAT(J3,1,3),ELLMAT(J3,2,1),ELLMAT(J3,2,2),ELLMAT(J3,2,3),
     &                       ELLMAT(J3,3,1),ELLMAT(J3,3,2),ELLMAT(J3,3,3),
     &            'ATOM_VECTOR',QMINP(J1,3*NATOMS/2+3*(J2-1)+1),
     &            QMINP(J1,3*NATOMS/2+3*(J2-1)+2),QMINP(J1,3*NATOMS/2+3*(J2-1)+3)
                END DO
               ELSE IF (GBT.OR.GBDT) THEN
                  CALL AATOEULER(QMINP(J1,3*NATOMS/2+3*(J2-1)+1),QMINP(J1,3*NATOMS/2+3*(J2-1)+2),
     &                           QMINP(J1,3*NATOMS/2+3*(J2-1)+3),EULERPHIDEG,EULERPSIDEG,EULERTHETADEG)

                  WRITE(MYUNIT2,'(A5,2X,3F20.10,2X,A8,6F15.8,2X,A11,3F15.8)') 'O',QMINP(J1,3*(J2-1)+1),
     &                  QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3),
     &            'ELLIPSE ',GBKAPPA,1.0D0,1.0D0,EULERPHIDEG,EULERPSIDEG,EULERTHETADEG,
     &            'ATOM_VECTOR',QMINP(J1,3*NATOMS/2+3*(J2-1)+1),
     &            QMINP(J1,3*NATOMS/2+3*(J2-1)+2),QMINP(J1,3*NATOMS/2+3*(J2-1)+3)
               ELSE IF (LJCAPSIDT) THEN
                  CALL AATOEULER(QMINP(J1,3*NATOMS/2+3*(J2-1)+1),QMINP(J1,3*NATOMS/2+3*(J2-1)+2),
     &                           QMINP(J1,3*NATOMS/2+3*(J2-1)+3),EULERPHIDEG,EULERPSIDEG,EULERTHETADEG)

                  WRITE(MYUNIT2,'(A5,2X,3F20.10,2X,A8,6F15.8,2X,A11,3F15.8)') 'O',QMINP(J1,3*(J2-1)+1),
     &                  QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3),
     &            'ELLIPSE ',1.0D0,1.0D0,1.0D0,EULERPHIDEG,EULERPSIDEG,EULERTHETADEG,
     &            'ATOM_VECTOR',QMINP(J1,3*NATOMS/2+3*(J2-1)+1),
     &            QMINP(J1,3*NATOMS/2+3*(J2-1)+2),QMINP(J1,3*NATOMS/2+3*(J2-1)+3)
               ENDIF
            ENDDO
         ENDDO
         CLOSE(MYUNIT2)

      ELSE IF (TIP) THEN
         OPEN(UNIT=26,FILE='TIP.XYZ',STATUS='UNKNOWN')
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
         OPEN(UNIT=26,FILE='CAPSID.XYZ',STATUS='UNKNOWN')
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

         OPEN(UNIT=26, FILE='DBP.XYZ', STATUS='UNKNOWN')
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
     &             'ATOM_VECTOR', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3) 
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

         OPEN(UNIT=26, FILE='GBE.XYZ', STATUS='UNKNOWN')
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

               WRITE(26,'(A5,2X,3F20.10,2X,A8,6F20.10)')
     &         'O', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3),
     &         'ELLIPSE', 2.D0*ESA(1), 2.D0*ESA(2), 2.D0*ESA(3), PHI, THT, 0.D0

            ENDDO

         ENDDO

         RETURN

      ELSE IF (LINRODT) THEN

         OPEN(UNIT=26, FILE='LINROD.XYZ', STATUS='UNKNOWN')
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
     &            'O', RBCOORDS(3*J4-2), RBCOORDS(3*J4-1), RBCOORDS(3*J4), 'ATOM_VECTOR', P(1), P(2), P(3)

               ENDDO

            ENDDO

         ENDDO

         CLOSE(UNIT=26)

         RETURN

      ELSE IF (LWOTPT) THEN

         OPEN(UNIT=26, FILE='LWOTP.XYZ', STATUS='UNKNOWN')
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
     &            'O', RBCOORDS(3*J4-2), RBCOORDS(3*J4-1), RBCOORDS(3*J4), 'ATOM_VECTOR', P(1), P(2), P(3)

               ENDDO

            ENDDO

         ENDDO

         CLOSE(UNIT=26)

         RETURN

      ELSE IF (MSTBINT) THEN

         OPEN(UNIT=26, FILE='MSTBIN.XYZ', STATUS='UNKNOWN')
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
     &               'O', RBCOORDS(3*J4-2), RBCOORDS(3*J4-1), RBCOORDS(3*J4), 'ATOM_VECTOR', P(1), P(2), P(3)
                  ELSE
                     WRITE(26,'(A4,3F20.10,2X,A12,2X,3F20.10)')
     &               'N', RBCOORDS(3*J4-2), RBCOORDS(3*J4-1), RBCOORDS(3*J4), 'ATOM_VECTOR', P(1), P(2), P(3)
                  ENDIF

               ENDDO

            ENDDO

         ENDDO

         CLOSE(UNIT=26)

         RETURN

      ELSE IF (MSSTOCKT) THEN

         OPEN(UNIT=26, FILE='MSSTOCK.XYZ', STATUS='UNKNOWN')
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
     &            'O', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3), 'ATOM_VECTOR', P(1), P(2), P(3)

               ENDDO

            ENDDO

         ENDDO

         CLOSE(UNIT=26)

         OPEN(UNIT=28, FILE='MSSTKTR.XYZ', STATUS='UNKNOWN')
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
     &            'O', RBCOORDS(3*J4-2), RBCOORDS(3*J4-1), RBCOORDS(3*J4), 'ATOM_VECTOR', P(1), P(2), P(3)

               ENDDO

            ENDDO

         ENDDO

         CLOSE(UNIT=28)

         RETURN

      ELSE IF (MULTPAHAT) THEN

         CALL VIEWMULTPAHA()
         RETURN

      ELSE IF (NPAHT .OR. PAHAT .OR. PAHW99T) THEN

         OPEN(UNIT=26, FILE='RIGID.XYZ', STATUS='UNKNOWN')
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

!|GD351>

      ELSE IF (PATCHY) THEN

         CALL VIEWPATCHY()
         RETURN

!<GD351|

      ELSE IF (PYGT .OR. PYGDPT) THEN

         OPEN(UNIT=26, FILE='ELLIPSOIDS.XYZ', STATUS='UNKNOWN')
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

               WRITE(26,'(A5,2X,3F20.10,2X,A8,6F20.10)')                 
     &              'O', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3),         
     &              'ELLIPSE', 2.D0*PYA1(1), 2.D0*PYA1(2), 2.D0*PYA1(3), PHI, THT, CHI 

            ENDDO

         ENDDO

         CLOSE(UNIT=26)

         RETURN

      ELSE IF (STOCKAAT) THEN

         OPEN(UNIT=26, FILE='STOCKAA.XYZ', STATUS='UNKNOWN')
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
     &         'ATOM_VECTOR', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)


            ENDDO

         ENDDO

         CLOSE (UNIT=26)

         RETURN

      ELSE IF (TDHDT) THEN

         CALL VIEWTDHD()
         RETURN

      ELSE IF (WATERDCT .OR. WATERKZT) THEN

         OPEN(UNIT=26, FILE='RIGID.XYZ', STATUS='UNKNOWN')
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
         OPEN(UNIT=26,FILE='RIGID.XYZ',STATUS='UNKNOWN')
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

      SUBROUTINE AMBERDUMP(J1,QMINP)
      USE COMMONS
      USE MODAMBER
      IMPLICIT NONE


      CHARACTER(LEN=25) COORDFILE
      CHARACTER(LEN=2) FNAME
      INTEGER J1
      DOUBLE PRECISION QMINP(NSAVE,3*NATOMS)

      IF (J1.LT.10) THEN
         WRITE (FNAME,'(I1)') J1
      ELSE
         WRITE (FNAME,'(I2)') J1
      ENDIF

      DO A=1,ATOMS
        X(A)=QMINP(J1,3*A-2)
        Y(A)=QMINP(J1,3*A-1)
        Z(A)=QMINP(J1,3*A)
      END DO

      COORDFILE='ACOORDS.DUMP.'//FNAME

      OPEN (UNIT=4,IOSTAT=IOS,FILE=COORDFILE,STATUS='UNKNOWN')

      DO A=1,ATOMS
        WRITE (UNIT=4,FMT='(A1,2X,A2,2X,I3,2X,I3,2X,F7.3,3X,F7.3,3X,F7.3)') LABEL(A),TYPECH(A),
     1        A,BONDEDTO(A),X(A),Y(A),Z(A)
      ENDDO

      WRITE (UNIT=4,FMT='(A3)') 'END'
      WRITE (UNIT=4,FMT='(A)') ' '
      WRITE (UNIT=4,FMT='(A4,7X,I2)') 'LOOP',RINGS

      DO A=1,RINGS
        WRITE (UNIT=4,FMT='(I3,4X,I3)') LOOPATOM(2*A-1),LOOPATOM(2*A)
      END DO

      WRITE (UNIT=4,FMT='(A)') ' '
      WRITE (UNIT=4,FMT='(A7)') 'CHARGES'

      DO A=1,ATOMS
        Q(A)=Q(A)/18.2223
        WRITE (UNIT=4,FMT='(I3,2X,F7.4)') A,Q(A)
      END DO

      WRITE (UNIT=4,FMT='(A3)') 'END'

      RETURN

      END
C
C  SUBROUTINE TO CONVERT CAPSID COFM AND DV COORDINATES TO PENATGONS.
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

      COORDS(1) =     C2A1 - C3A1*L12 + X1
      COORDS(2) =     -(C3A1*L1*M1) - N1*S1 + Y1
      COORDS(3) =     -(C3A1*L1*N1) + M1*S1 + Z1
      COORDS(4) =     C2A1*NUM4 - C3A1*L1*(M1*NUM3 + L1*NUM4) + N1*NUM3*S1 + X1
      COORDS(5) =     C2A1*NUM3 - C3A1*M1*(M1*NUM3 + L1*NUM4) - N1*NUM4*S1 + Y1
      COORDS(6) =     -(C3A1*N1*(M1*NUM3 + L1*NUM4)) - L1*NUM3*S1 + M1*NUM4*S1 + Z1
      COORDS(7) =     C2A1*NUM1 - C3A1*L1*(L1*NUM1 + M1*NUM2) + N1*NUM2*S1 + X1
      COORDS(8) = C2A1*NUM2 - C3A1*M1*(L1*NUM1 + M1*NUM2) - N1*NUM5*S1 + Y1
      COORDS(9) = -(C3A1*N1*(L1*NUM1 + M1*NUM2)) + M1*NUM1*S1 - L1*NUM2*S1 + Z1
      COORDS(10) = C2A1*NUM1 + C3A1*L1*(-(L1*NUM1) + M1*NUM2) - N1*NUM2*S1 + X1
      COORDS(11) = -(C2A1*NUM2) + C3A1*M1*(-(L1*NUM1) + M1*NUM2) - N1*NUM5*S1 + Y1
      COORDS(12) = -(C3A1*L1*N1*NUM1) + C3A1*M1*N1*NUM2 + M1*NUM1*S1 + L1*NUM2*S1 + Z1
      COORDS(13) = C2A1*NUM4 + C3A1*L1*(M1*NUM3 - L1*NUM4) - N1*NUM3*S1 + X1
      COORDS(14) = -(C2A1*NUM3) + C3A1*M1*(M1*NUM3 - L1*NUM4) - N1*NUM4*S1 + Y1
      COORDS(15) = C3A1*N1*(M1*NUM3 - L1*NUM4) + L1*NUM3*S1 + M1*NUM4*S1 + Z1
C     COORDS(16)= (-(C3A1*L1*N1) - M1*S1 + 2*X1)/2.
C     COORDS(17)= -(C3A1*M1*N1)/2. + (L1*S1)/2. + Y1
C     COORDS(18)= (C2A1 - C3A1*N12 + 2*Z1)/2.
      COORDS(16)= -(C3A1*HEIGHT*L1*N1) - HEIGHT*M1*S1 + X1
      COORDS(17)= -(C3A1*HEIGHT*M1*N1) + HEIGHT*L1*S1 + Y1
      COORDS(18)= C2A1*HEIGHT - C3A1*HEIGHT*N12 + Z1

      RETURN
      END
C
C  SUBROUTINE TO CONVERT RIGID BODY COFM AND DV COORDINATES TO MOLECULAR SITES.
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
         COORDS(3*(J1-1)+1)=C2A1*SITE(J1,1) + S1*(N1*SITE(J1,2) - M1*SITE(J1,3)) - 
     1                      C3A1*L1*(L1*SITE(J1,1) + M1*SITE(J1,2) + N1*SITE(J1,3)) + X1
         COORDS(3*(J1-1)+2)=C2A1*SITE(J1,2) + S1*(-(N1*SITE(J1,1)) + L1*SITE(J1,3)) 
     1                    - C3A1*M1*(L1*SITE(J1,1) + M1*SITE(J1,2) + N1*SITE(J1,3)) + Y1
         COORDS(3*(J1-1)+3)=S1*(M1*SITE(J1,1) - L1*SITE(J1,2)) + C2A1*SITE(J1,3) 
     1                    - C3A1*N1*(L1*SITE(J1,1) + M1*SITE(J1,2) + N1*SITE(J1,3)) + Z1
      ENDDO 

      RETURN
      END
C
C  SUBROUTINE TO CONVERT TIP OXYGEN AND DV COORDINATES TO CARTESIANS.
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
      COORDS(4) = 0.756950327*C2A1 - C3A1*L1*(0.756950327*L1 - 0.585882276*N1) + 0.585882276*M1*S1 + X1
      COORDS(5) = -(C3A1*M1*(0.756950327*L1 - 0.585882276*N1)) + (-0.585882276*L1 - 0.756950327*N1)*S1 + Y1
      COORDS(6) = -0.585882276*C2A1 - C3A1*(0.756950327*L1 - 0.585882276*N1)*N1 + 0.756950327*M1*S1 + Z1
      COORDS(7) = -0.756950327*C2A1 + C3A1*L1*(0.756950327*L1 + 0.585882276*N1) + 0.585882276*M1*S1 + X1
      COORDS(8) = C3A1*M1*(0.756950327*L1 + 0.585882276*N1) + (-0.585882276*L1 + 0.756950327*N1)*S1 + Y1
      COORDS(9) = -0.585882276*C2A1 + C3A1*(0.756950327*L1 + 0.585882276*N1)*N1 - 0.756950327*M1*S1 + Z1

      RETURN
      END

