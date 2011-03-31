C   OPTIM: A program for optimizing geometries and calculating reaction pathways
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of OPTIM.
C
C   OPTIM is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   OPTIM is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
      PROGRAM OPTIM4

      USE COMMONS
      USE PORFUNCS
      !USE OPTIMHEADER
      USE KEY, ONLY: FILTHSTR,SEQ,NUMGLY,TARFL,CASTEPJOB,CP2KJOB,ONETEPJOB
      USE MODAMBER9, ONLY: AMBERSTR,AMBERSTR1,INPCRD,ATMASS1

      IMPLICIT NONE

      INTEGER ITEM, NITEMS, LOC, LINE, NCR, NERROR, IR, LAST, J1
      LOGICAL VARIABLES,CASTEP,ONETEP,CP2K,DFTP,CPMD,END,CAT,SKIPBL,CLEAR,ECHO,AMBER,AMBERT,NABT,RINGPOLYMERT
      COMMON /BUFINF/ ITEM, NITEMS, LOC(80), LINE, SKIPBL, CLEAR, NCR,
     &                NERROR, IR, ECHO, LAST, CAT
      DOUBLE PRECISION DUMMY
      CHARACTER ZDUM*5
      INTEGER J2, NELEMENTS, LSYS, NTYPE(105), IOS, NARG, FILTH, FILTH2
      CHARACTER FNAME*80, TSTRING*80
      CHARACTER(LEN=80) :: SYS
      CHARACTER WORD*16
      CHARACTER(LEN=10)  check
      CHARACTER(LEN=20) OTEMP, OSTRING, CSTRING
      CHARACTER(LEN=21) DSTRING1, DSTRING2
      CHARACTER(LEN=80) ARGSTRING, MYLINE
      LOGICAL AMH
      INTEGER :: NRES,I_RES,NOGLY

      !CALL PRINTHEADER
      CASTEP=.FALSE.
      ONETEP=.FALSE.
      CP2K=.FALSE. 
      CPMD=.FALSE.
      VARIABLES=.FALSE.
      RINGPOLYMERT=.FALSE.
      AMBER=.FALSE.
      AMBERT=.FALSE.
      NABT=.FALSE.

      CALL READ_CMD_ARGS

      IF (FILTH2.EQ.0) THEN
         OPEN(5,FILE='odata',STATUS='OLD')
      ELSE
         WRITE(OTEMP,*) FILTH2
         WRITE(OSTRING,'(A)') 'odata.' // TRIM(ADJUSTL(OTEMP))
         OPEN(5,FILE=OSTRING,STATUS='OLD')
      ENDIF

190   CALL INPUT(END)
      IF (.NOT. END) CALL READU(WORD)
      IF (END.OR.WORD.EQ.'STOP'.OR.WORD.EQ.'POINTS') GOTO 200
      IF (WORD.EQ.'    ' .OR.WORD.EQ.'NOTE'.OR.WORD.EQ.'COMMENT'.OR. WORD .EQ. '\\') THEN
         GOTO 190
      ELSE IF ((WORD.EQ.'CPMD').OR.(WORD.EQ.'CPMDC')) THEN
         CPMD=.TRUE.
         IF (NITEMS.GT.1) THEN
            CALL READA(SYS)
         ELSE
            WRITE(*,'(A)') ' ERROR - no CPMD system specified'
            STOP
         ENDIF
         DO J1=1,80
            IF (SYS(J1:J1).EQ.' ') THEN
               LSYS=J1-1
               GOTO 112
            ENDIF
         ENDDO
112      CONTINUE
      ELSE IF ((WORD.EQ.'ONETEP').OR.(WORD.EQ.'ONETEPC')) THEN
         ONETEP=.TRUE. 
         IF (WORD.EQ.'ONETEP') DFTP=.TRUE.
         IF (NITEMS.GT.2) THEN
            CALL READA(ONETEPJOB)
            CALL READA(SYS)
         ELSE
            WRITE(*,'(A)') 'getparams> ERROR - ONETEP input mangled'
            STOP
         ENDIF
         DO J1=1,80
            IF (SYS(J1:J1).EQ.' ') THEN
               LSYS=J1-1
               GOTO 231
            ENDIF
         ENDDO
231      CONTINUE
      ELSE IF ((WORD.EQ.'CP2K').OR.(WORD.EQ.'CP2KC')) THEN 
         CP2K=.TRUE.
         IF (WORD.EQ.'CP2K') DFTP=.TRUE. 
         IF (NITEMS.GT.2) THEN
            CALL READA(CP2KJOB)
            CALL READA(SYS)
         ELSE
            WRITE(*,'(A)') 'getparams> ERROR - CP2K input mangled'
            STOP
         ENDIF
         DO J1=1,80
            IF (SYS(J1:J1).EQ.' ') THEN
               LSYS=J1-1
               GOTO 271
            ENDIF
         ENDDO
271      CONTINUE
      ELSE IF ((WORD.EQ.'CASTEP').OR.(WORD.EQ.'CASTEPC')) THEN
         CASTEP=.TRUE.
         IF (WORD.EQ.'CASTEP') DFTP=.TRUE.
         IF (NITEMS.GT.2) THEN
            CALL READA(CASTEPJOB)
            CALL READA(SYS)
         ELSE
            WRITE(*,'(A)') 'getparams> ERROR - CASTEP input mangled'
            STOP
         ENDIF
         DO J1=1,80
            IF (SYS(J1:J1).EQ.' ') THEN
               LSYS=J1-1
               GOTO 211
            ENDIF
         ENDDO
211      CONTINUE
      ELSE IF (WORD.EQ.'RINGPOLYMER') THEN
         RINGPOLYMERT=.TRUE.
         GOTO 200
      ELSE IF (WORD.EQ.'VARIABLES') THEN
         VARIABLES=.TRUE.
         GOTO 200
      ELSE IF (WORD.EQ.'AMBER') THEN
         AMBER=.TRUE.
         NATOMS=0
         OPEN (UNIT=9,FILE='coords.amber',STATUS='OLD')
310      READ (UNIT=9,IOSTAT=ios,FMT='(A3)') check
         IF (ios.LT.0) THEN
           PRINT *,'End of file before all information specified'
           STOP
         ENDIF
         IF (check.EQ.'end' .OR. check.EQ.'END' .OR. check.EQ.'End') THEN
            CLOSE(9)
            GOTO 400
         ENDIF
         NATOMS=NATOMS+1
         GOTO 310
      ELSE IF (WORD.EQ.'AMBER9'.OR.(WORD.EQ.'NAB')) THEN
         IF(WORD.EQ.'AMBER9') AMBERT=.TRUE.
         IF(WORD.EQ.'NAB') NABT=.TRUE.
!         NATOMS=0
        inpcrd='coords.inpcrd'
       IF(NITEMS<3) then
         CALL READA(amberstr)
          IF (FILTH2.NE.0) THEN
            WRITE(OTEMP,*) FILTH2
            write(ostring,'(A)') trim(adjustl(amberstr))//'.'//trim(adjustl(otemp))
            WRITE(*,*) 'ostring=', ostring
          ELSE
            write(ostring,'(A)') trim(adjustl(amberstr))
          END IF
          WRITE(*,'(A)') ' getparams> input coordinates for AMBER9 system will be read from ',trim(adjustl(amberstr)),ostring
         call amberinterface(natoms,2,inpcrd,6)
         CALL amber_readcoords(ostring)
       ELSE IF(NITEMS==3) then
         CALL READA(amberstr)
         CALL READA(amberstr1)
          WRITE(*,'(A)') ' getparams> input coordinates for AMBER9 system will be read from ', trim(adjustl(amberstr)),
     &                         'type: ', trim(adjustl(amberstr1))
          IF(trim(adjustl(amberstr1)).EQ.'inpcrd') then
           inpcrd=amberstr
           WRITE(*,'(A)') ' getparams> reading AMBER inpcrd coordinate format'
          ELSE
           WRITE(*,'(A)') ' getparams> ERROR - no other types defined currently than inpcrd'
          END IF
           call amberinterface(natoms,2,inpcrd,6)
       END IF
        GOTO 400           ! finished counting atoms, go to the end of the subroutine
      ELSE IF (WORD.EQ.'AMH') THEN

         AMH=.TRUE.
         NATOMS=0
         WRITE(6,*)'Entering GETPARAMS'

         OPEN(UNIT=30,FILE='pro.list',STATUS='OLD',FORM='FORMATTED')
         READ (30,1000)TARFL
1000     FORMAT(A5)
         CLOSE(30)

         OPEN(30,FILE='proteins/'//TARFL,STATUS='OLD')
            READ(30,*)
            READ(30,*)NRES
            IF (NRES.GT.500) THEN
                WRITE(6,*) 'FAILURE NRES GR THAN 500 COUNTATOMS'
                STOP
            ENDIF
            READ (30,25)(SEQ(I_RES),I_RES=1,NRES)
25         FORMAT(25(I2,1X))
          CLOSE(30)
               
          NOGLY = 0
          NUMGLY = 0
           DO I_RES=1,NRES
             IF (SEQ(I_RES).NE.8) NOGLY = NOGLY +1
             IF (SEQ(I_RES).EQ.8) NUMGLY = NUMGLY +1
           ENDDO
            NATOMS = NOGLY*3 + NUMGLY*2
            WRITE(6,*)'NATOMS NOGLY  NUMGLY  ',NATOMS,NOGLY,NUMGLY
           GOTO 400

       ELSE IF (WORD.EQ.'CHARMM') THEN
! DAE We are going to assume that there is a charmm format file in the directory called input.crd, and the first
! line of it will tell us the number of atoms. In 99% of old (OPTIM<3) runs this is what I did, but the old versions
! were more flexible in that any filename could be specified in the CHARMM bit of the odata file

         IF (FILTH2.NE.0) THEN
            WRITE(OTEMP,*) FILTH2
            WRITE(OSTRING,'(A)') 'input.crd.' // TRIM(ADJUSTL(OTEMP))
         ELSE
            WRITE(OSTRING,'(A)') 'input.crd'
         ENDIF
         OPEN (UNIT=9,FILE=OSTRING,STATUS='OLD',IOSTAT=ios)

         IF (ios /= 0) THEN
            WRITE(OSTRING,'(A)') 'input.crd'
            OPEN (UNIT=9,FILE=OSTRING,STATUS='OLD',IOSTAT=ios)
            if (ios == 0) THEN
         else
            WRITE(*,'(2A)') 'Thanks to our new dynamic memory allocation overlords, there must be a charmm-format file called ',
     &    '"input.crd" for CHARMM to find out the number of atoms. Feel free to recode to enable any filename to work properly.'
            STOP
         ENDIF
         ENDIF
         do
              read(9,*) myline
              if (myline(1:1)=='*') then ! SAT This is the goddamn CHARMM comment line
                    cycle
              else
                    read(myline,*) NATOMS
                    exit
              endif
         enddo
         CLOSE(9)

! DAE We also need to find out what MAXAIM is in CHARMM, and set MXATMS in OPTIM to be the same, so that those arrays which
! are passed between the two can be declared correctly. MXATMS is now stored in modmxatms.
       
         CALL GETMAXAIM

         GOTO 400

      ELSE IF (WORD.EQ.'UNRES') THEN
! jmc We are going to assume that there is a coords file. The first line of it will tell us the number of atoms.

         IF (FILTH2.EQ.0) THEN
            OPEN (UNIT=9,FILE='coords',STATUS='OLD',IOSTAT=ios)
         ELSE
            WRITE(CSTRING,'(A)') 'coords.'//TRIM(ADJUSTL(OTEMP))
            OPEN (UNIT=9,FILE=CSTRING,STATUS='OLD',IOSTAT=ios)
         ENDIF
         IF (ios /= 0) THEN
            WRITE(*,'(2A)') 'Thanks to our new dynamic memory allocation overlords, there must be a coords file present ',
     &    ' for OPTIM3 to find out the number of atoms for UNRES. Please make it so!'
            STOP
         ENDIF
         READ(9,*) NATOMS
         CLOSE(9)

         GOTO 400

      ENDIF
      GOTO 190

200   CONTINUE

      IF (CASTEP) THEN
         FNAME=SYS(1:LSYS) // '.cell'
         WRITE(*,'(A,A)') ' getparams> Counting atoms in file ',FNAME
         OPEN(UNIT=7,FILE=FNAME,STATUS='OLD')
         NATOMS=0
         coordsloop: DO
            READ(7,'(A21)') DSTRING1
            CALL UPPERCASE(DSTRING1)
C           WRITE(*,'(A21)') DSTRING1
            IF (DSTRING1(1:16).EQ.'%BLOCK POSITIONS') THEN
               DO 
                  READ(7,'(A21)') DSTRING1
                  CALL UPPERCASE(DSTRING1)
C                 WRITE(*,*) DSTRING1
                  IF (DSTRING1(1:19).EQ.'%ENDBLOCK POSITIONS') EXIT coordsloop
                  NATOMS=NATOMS+1
               ENDDO
            ENDIF
         ENDDO coordsloop
C        WRITE(*,'(A,I4,A)') 'getparams> CASTEP run for ',NATOMS,' atoms'
         CLOSE(7)
      ELSEIF (ONETEP) THEN
         FNAME=SYS(1:LSYS) // '.dat'
         WRITE(*,'(A,A)') 'getparams> Counting atoms in file ',FNAME
         OPEN(UNIT=7,FILE=FNAME,STATUS='OLD')
         NATOMS=0
         coordsloop2: DO
            READ(7,'(A21)') DSTRING1
            CALL UPPERCASE(DSTRING1)
C           WRITE(*,'(A21)') DSTRING1
            IF (DSTRING1(1:16).EQ.'%BLOCK POSITIONS') THEN
               DO 
                  READ(7,'(A21)') DSTRING1
                  CALL UPPERCASE(DSTRING1)
C                 WRITE(*,*) DSTRING1
                  IF (DSTRING1(1:19).EQ.'%ENDBLOCK POSITIONS') EXIT coordsloop2
                  NATOMS=NATOMS+1
               ENDDO
            ENDIF
         ENDDO coordsloop2
         WRITE(*,'(A,I4,A)') 'getparams> ONETEP run for ',NATOMS,' atoms'
         CLOSE(7)
      ELSEIF (CP2K) THEN 
         FNAME=SYS(1:LSYS) // '.inp'
         WRITE(*,'(A,A)') ' getparams> Counting atoms in file ',FNAME
         OPEN(UNIT=7,FILE=FNAME,STATUS='OLD')
         NATOMS=0
         coordsloop3: DO
            READ(7,'(A21)') DSTRING1
            DSTRING1=ADJUSTL(DSTRING1)
            CALL UPPERCASE(DSTRING1)
C           WRITE(*,'(A21)') DSTRING1
            IF (DSTRING1(1:6).EQ.'&COORD') THEN
               DO 
                  READ(7,'(A21)') DSTRING1
                  DSTRING1=ADJUSTL(DSTRING1)
                  CALL UPPERCASE(DSTRING1)
C                 WRITE(*,*) DSTRING1
                  IF (DSTRING1(1:10).EQ.'&END COORD') EXIT coordsloop3
                  IF (DSTRING1(1:6).EQ.'SCALED') NATOMS=NATOMS-1
                  IF (DSTRING1(1:4).EQ.'UNIT') NATOMS=NATOMS-1
                  NATOMS=NATOMS+1
               ENDDO
            ENDIF
         ENDDO coordsloop3
         WRITE(*,'(A,I4,A)') ' getparams> CP2K run for ',NATOMS,' atoms'
         CLOSE(7) 
      ELSE IF (CPMD) THEN
         FNAME=SYS(1:LSYS)
         WRITE(*,'(A,A)') 'getparams> Reading coordinates from file ',FNAME
         OPEN(UNIT=7,FILE=FNAME,STATUS='OLD')
         NATOMS=0
11       READ(7,'(A)') FNAME
         IF (FNAME(1:6).EQ.'&ATOMS') THEN
            J1=0
12          READ(7,'(A)') TSTRING
            IF (TSTRING(1:1).EQ.'*') THEN
               J1=J1+1
               READ(7,'(A)') FNAME
               READ(7,*) NTYPE(J1)
               DO J2=1,NTYPE(J1)
                  READ(7,*) DUMMY
               ENDDO
               NATOMS=NATOMS+NTYPE(J1)
               GOTO 12
            ELSE IF (TSTRING(1:1).EQ.' ') THEN
               GOTO 12
            ELSE IF (TSTRING(1:4).EQ.'&END') THEN
               GOTO 13
            ENDIF
         ELSE
            GOTO 11
         ENDIF
13       CONTINUE
      ELSE IF (VARIABLES) THEN
         NATOMS=0
110      CALL INPUT(END)
         IF (.NOT.END) THEN
            CALL READF(DUMMY)
            NATOMS=NATOMS+1
            GOTO 110
         ENDIF
      ELSE IF (RINGPOLYMERT) THEN
         NATOMS=0
111      CALL INPUT(END)
         IF (.NOT.END) THEN
            CALL READF(DUMMY)
            NATOMS=NATOMS+NITEMS
            GOTO 111
         ENDIF
      ELSE 
         NATOMS=0
300      CALL INPUT(END)
         IF (.NOT.END) THEN
            CALL READA(ZDUM)
            IF (ZDUM(1:2).EQ.'  ') GOTO 300
            CALL READF(DUMMY)
            CALL READF(DUMMY)
            CALL READF(DUMMY)
            NATOMS=NATOMS+1
            GOTO 300
         ENDIF
      ENDIF
400   CLOSE(5)

      WRITE(*,'(A,I6)') ' getparams> Number of atoms (or variables)  determined as ',NATOMS

      CALL OPTIM(FILTH,FILTH2,ARGSTRING)

      STOP
      END
